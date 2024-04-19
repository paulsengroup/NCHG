// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public License along
// with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#pragma once

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/expected_values_aggregator.hpp>
#include <hictk/pixel.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <hictk/transformers/pixel_merger.hpp>
#include <highfive/H5File.hpp>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/mad_max_filter.hpp"

namespace nchg {

template <typename File>
inline ExpectedValues<File>::ExpectedValues(std::shared_ptr<const File> file, double mad_max,
                                            std::uint64_t min_delta, std::uint64_t max_delta)
    : _fp(std::move(file)), _mad_max(mad_max), _min_delta(min_delta), _max_delta(max_delta) {
  if (_mad_max < 0 || !std::isfinite(_mad_max)) {
    throw std::logic_error("mad_max_ should be a non-negative value");
  }
  if (_min_delta >= _max_delta) {
    throw std::logic_error("min_delta should be strictly less than max_delta");
  }
  if (_fp) {
    compute_expected_values_cis();
    compute_expected_values_trans();
  }
}

template <typename File>
inline ExpectedValues<File> ExpectedValues<File>::cis_only(std::shared_ptr<const File> file,
                                                           double mad_max, std::uint64_t min_delta,
                                                           std::uint64_t max_delta) {
  ExpectedValues ev(nullptr, mad_max, min_delta, max_delta);
  ev._fp = std::move(file);
  if (ev._fp) {
    ev.compute_expected_values_cis();
  }
  return ev;
}

template <typename File>
inline ExpectedValues<File> ExpectedValues<File>::trans_only(std::shared_ptr<const File> file,
                                                             double mad_max) {
  ExpectedValues ev(nullptr, mad_max);
  ev._fp = std::move(file);
  if (ev._fp) {
    ev.compute_expected_values_trans();
  }
  return ev;
}

template <typename File>
inline ExpectedValues<File> ExpectedValues<File>::chromosome_pair(std::shared_ptr<const File> file,
                                                                  const hictk::Chromosome &chrom1,
                                                                  const hictk::Chromosome &chrom2,
                                                                  double mad_max,
                                                                  std::uint64_t min_delta,
                                                                  std::uint64_t max_delta) {
  SPDLOG_INFO(FMT_STRING("computing expected values for {}:{}..."), chrom1.name(), chrom2.name());

  ExpectedValues ev(nullptr, mad_max, min_delta, max_delta);
  ev._fp = std::move(file);
  if (chrom1 == chrom2) {
    ev.compute_expected_values_cis();
    return ev;
  }

  const auto sel = ev._fp->fetch(chrom1.name(), chrom2.name());
  const hictk::transformers::JoinGenomicCoords jsel(
      sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), ev._fp->bins_ptr());

  if (!ev._bin_masks.contains(std::make_pair(chrom1, chrom2))) {
    ev.add_bin_mask(chrom1, chrom2,
                    mad_max_filtering(jsel.begin(), jsel.end(), chrom1, chrom2,
                                      ev._fp->resolution(), ev._mad_max));
  }

  const ExpectedMatrix em(jsel.begin(), jsel.end(), chrom1, chrom2, ev._fp->bins(),
                          std::vector<double>{}, 0, *ev.bin_mask(chrom1, chrom2).first,
                          *ev.bin_mask(chrom1, chrom2).second,
                          std::numeric_limits<std::uint64_t>::max());

  ev._expected_values_trans.emplace(std::make_pair(chrom1, chrom2), em.nnz_avg());

  return ev;
}

template <typename File>
inline ExpectedValues<File> ExpectedValues<File>::deserialize(const std::filesystem::path &path) {
  ExpectedValues<File> ev{nullptr};
  HighFive::File f(path.string());

  auto [weights, scaling_factors] = deserialize_cis_profiles(f);
  ev._expected_weights = std::move(weights);
  ev._expected_scaling_factors = std::move(scaling_factors);
  ev._expected_values_trans = deserialize_trans_profiles(f);

  return ev;
}

template <typename File>
inline const std::vector<double> &ExpectedValues<File>::weights() const noexcept {
  return _expected_weights;
}

template <typename File>
inline double ExpectedValues<File>::mad_max() const noexcept {
  return _mad_max;
}
template <typename File>
inline std::uint64_t ExpectedValues<File>::min_delta() const noexcept {
  return _min_delta;
}
template <typename File>
inline std::uint64_t ExpectedValues<File>::max_delta() const noexcept {
  return _max_delta;
}

template <typename File>
inline std::shared_ptr<const std::vector<bool>> ExpectedValues<File>::bin_mask(
    const hictk::Chromosome &chrom) const {
  return _bin_masks.at(std::make_pair(chrom, chrom)).first;
}

template <typename File>
inline std::pair<std::shared_ptr<const std::vector<bool>>, std::shared_ptr<const std::vector<bool>>>
ExpectedValues<File>::bin_mask(const hictk::Chromosome &chrom1,
                               const hictk::Chromosome &chrom2) const {
  return _bin_masks.at(std::make_pair(chrom1, chrom2));
}

template <typename File>
inline std::vector<double> ExpectedValues<File>::expected_values(const hictk::Chromosome &chrom,
                                                                 bool rescale) const {
  if (_expected_weights.empty()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("expected values for \"{}\" are not available: out of range"), chrom.name()));
  }

  const auto num_bins = (chrom.size() + _fp->resolution() - 1) / _fp->resolution();
  assert(num_bins <= _expected_weights.size());

  std::vector<double> weights(_expected_weights.begin(),
                              _expected_weights.begin() + static_cast<std::ptrdiff_t>(num_bins));

  if (!rescale) {
    return weights;
  }

  const auto sf = _expected_scaling_factors.at(chrom);
  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [&](const auto n) { return n / sf; });

  return weights;
}

template <typename File>
inline double ExpectedValues<File>::expected_value(const hictk::Chromosome &chrom1,
                                                   const hictk::Chromosome &chrom2) const {
  return _expected_values_trans.at(EVTKey{chrom1, chrom2});
}

template <typename File>
inline const phmap::btree_map<hictk::Chromosome, double> &ExpectedValues<File>::scaling_factors()
    const noexcept {
  return _expected_scaling_factors;
}

template <typename File>
inline double ExpectedValues<File>::scaling_factor(const hictk::Chromosome &chrom) const {
  return scaling_factors().at(chrom);
}

template <typename File>
inline auto ExpectedValues<File>::expected_matrix(const hictk::Chromosome &chrom,
                                                  const hictk::BinTable &bins, PixelIt first_pixel,
                                                  PixelIt last_pixel) const
    -> ExpectedMatrix<PixelIt> {
  return {std::move(first_pixel),
          std::move(last_pixel),
          chrom,
          chrom,
          bins,
          _expected_weights,
          _expected_scaling_factors.at(chrom),
          *bin_mask(chrom),
          *bin_mask(chrom),
          _min_delta,
          _max_delta};
}

template <typename File>
inline auto ExpectedValues<File>::expected_matrix(const hictk::Chromosome &chrom) const
    -> ExpectedMatrix<PixelIt> {
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }
  const auto sel = _fp->fetch(chrom.name());
  const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<N>(), sel.template end<N>(),
                                                    _fp->bins_ptr());
  return expected_matrix(chrom, _fp->bins(), jsel.begin(), jsel.end());
}

template <typename File>
inline auto ExpectedValues<File>::expected_matrix(const hictk::Chromosome &chrom1,
                                                  const hictk::Chromosome &chrom2,
                                                  const hictk::BinTable &bins, PixelIt first_pixel,
                                                  PixelIt last_pixel) const
    -> ExpectedMatrix<PixelIt> {
  if (chrom1 == chrom2) {
    return expected_matrix(chrom1, bins, std::move(first_pixel), std::move(last_pixel));
  }
  return {std::move(first_pixel),
          std::move(last_pixel),
          chrom1,
          chrom2,
          bins,
          std::vector<double>{},
          1.0,
          *bin_mask(chrom1, chrom2).first,
          *bin_mask(chrom1, chrom2).second,
          _min_delta,
          _max_delta};
}

template <typename File>
inline auto ExpectedValues<File>::expected_matrix(const hictk::Chromosome &chrom1,
                                                  const hictk::Chromosome &chrom2) const
    -> ExpectedMatrix<PixelIt> {
  if (chrom1 == chrom2) {
    return expected_matrix(chrom1);
  }

  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }
  const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
  const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<N>(), sel.template end<N>(),
                                                    _fp->bins_ptr());
  return expected_matrix(chrom1, chrom2, _fp->bins(), jsel.begin(), jsel.end());
}

template <typename File>
inline void ExpectedValues<File>::serialize(const std::filesystem::path &path) const {
  SPDLOG_INFO(FMT_STRING("writing expected value profiles to {}..."), path);
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  HighFive::File f(path.string(), HighFive::File::Create);

  const auto source_file = std::filesystem::path{_fp->path()}.filename().string();
  f.createAttribute("source-file", source_file);
  f.createAttribute("resolution", _fp->resolution());

  serialize_chromosomes(f, _fp->chromosomes());
  serialize_cis_profiles(f, _expected_weights, _expected_scaling_factors);
  serialize_trans_profiles(f, _expected_values_trans);
}

template <typename File>
inline void ExpectedValues<File>::compute_expected_values_cis() {
  SPDLOG_INFO(FMT_STRING("initializing expected matrix weights from genome-wide interactions..."));
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  using PixelSel = decltype(std::declval<File>().fetch("chr1"));

  std::vector<ThinPixelIt> heads{};
  std::vector<ThinPixelIt> tails{};

  std::vector<PixelSel> sels{};

  for (const auto &chrom : _fp->chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    auto sel = _fp->fetch(chrom.name());

    auto first = sel.template begin<N>();
    auto last = sel.template end<N>();

    if (first != last) {
      if (!_bin_masks.contains(std::make_pair(chrom, chrom))) {
        const hictk::transformers::JoinGenomicCoords jsel(first, last, _fp->bins_ptr());
        add_bin_mask(
            chrom, mad_max_filtering(jsel.begin(), jsel.end(), chrom, _fp->resolution(), _mad_max));
      }

      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
      sels.emplace_back(std::move(sel));
    }
  }

  if (heads.empty()) {
    _expected_weights.resize(_fp->bins().size(), 0);
    for (const auto &chrom : _fp->chromosomes()) {
      _expected_scaling_factors.emplace(chrom, 1.0);
    }
    return;
  }

  hictk::transformers::PixelMerger merger(std::move(heads), std::move(tails));
  const hictk::transformers::JoinGenomicCoords mjsel(merger.begin(), merger.end(), _fp->bins_ptr());

  hictk::ExpectedValuesAggregator aggr(_fp->bins_ptr());
  std::for_each(mjsel.begin(), mjsel.end(), [&](const hictk::Pixel<N> &p) {
    const auto &mask = *bin_mask(p.coords.bin1.chrom());
    const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
    const auto bin1_id = p.coords.bin1.rel_id();
    const auto bin2_id = p.coords.bin2.rel_id();
    const auto bin1_masked = !mask.empty() && mask[bin1_id];
    const auto bin2_masked = !mask.empty() && mask[bin2_id];
    if (delta >= _min_delta && delta < _max_delta && !bin1_masked && !bin2_masked) {
      aggr.add(p);
    }
  });
  aggr.compute_density();

  _expected_scaling_factors = aggr.scaling_factors();
  for (const auto &chrom : _fp->chromosomes()) {
    if (!_expected_scaling_factors.contains(chrom)) {
      _expected_scaling_factors.emplace(chrom, std::numeric_limits<double>::quiet_NaN());
    }
  }

  _expected_weights = aggr.weights();
  const auto &chrom = _fp->chromosomes().longest_chromosome();
  const auto num_bins = (chrom.size() + _fp->resolution() - 1) / _fp->resolution();
  _expected_weights.resize(num_bins, std::numeric_limits<double>::quiet_NaN());
}

template <typename File>
inline void ExpectedValues<File>::compute_expected_values_trans() {
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }
  const auto &chroms = _fp->chromosomes();
  for (std::uint32_t chrom1_id = 0; chrom1_id < chroms.size(); ++chrom1_id) {
    const auto &chrom1 = chroms.at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < chroms.size(); ++chrom2_id) {
      const auto &chrom2 = chroms.at(chrom2_id);

      SPDLOG_INFO(FMT_STRING("processing {}:{}..."), chrom1.name(), chrom2.name());

      const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
      const hictk::transformers::JoinGenomicCoords jsel(
          sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), _fp->bins_ptr());

      if (!_bin_masks.contains(std::make_pair(chrom1, chrom2))) {
        add_bin_mask(chrom1, chrom2,
                     mad_max_filtering(jsel.begin(), jsel.end(), chrom1, chrom2, _fp->resolution(),
                                       _mad_max));
      }

      const ExpectedMatrix em(jsel.begin(), jsel.end(), chrom1, chrom2, _fp->bins(),
                              std::vector<double>{}, 0, *bin_mask(chrom1, chrom2).first,
                              *bin_mask(chrom1, chrom2).second,
                              std::numeric_limits<std::uint64_t>::max());
      _expected_values_trans.emplace(std::make_pair(chrom1, chrom2), em.nnz_avg());
    }
  }
}

template <typename File>
inline void ExpectedValues<File>::add_bin_mask(const hictk::Chromosome &chrom,
                                               std::vector<bool> &&mask) {
  const auto key = std::make_pair(chrom, chrom);
  if (_bin_masks.contains(key)) {
    return;
  }
  auto value = std::make_shared<const std::vector<bool>>(std::move(mask));
  _bin_masks.emplace(key, std::make_pair(value, value));
}

template <typename File>
void ExpectedValues<File>::add_bin_mask(const hictk::Chromosome &chrom1,
                                        const hictk::Chromosome &chrom2,
                                        std::pair<std::vector<bool>, std::vector<bool>> &&masks) {
  const auto key = std::make_pair(chrom1, chrom2);
  if (_bin_masks.contains(key)) {
    return;
  }
  _bin_masks.emplace(
      key, std::make_pair(std::make_shared<const std::vector<bool>>(std::move(masks.first)),
                          std::make_shared<const std::vector<bool>>(std::move(masks.second))));
}

template <typename File>
inline void ExpectedValues<File>::serialize_chromosomes(HighFive::File &f,
                                                        const hictk::Reference &chroms) {
  auto grp = f.createGroup("chroms");

  std::vector<std::string> chromosome_names(chroms.size());
  std::transform(chroms.begin(), chroms.end(), chromosome_names.begin(),
                 [](const hictk::Chromosome &c) { return std::string{c.name()}; });
  grp.createDataSet("name", chromosome_names);

  std::vector<std::uint32_t> chromosome_sizes(chroms.size());
  std::transform(chroms.begin(), chroms.end(), chromosome_sizes.begin(),
                 [](const hictk::Chromosome &c) { return c.size(); });
  grp.createDataSet("length", chromosome_sizes);
}

template <typename File>
inline void ExpectedValues<File>::serialize_cis_profiles(
    HighFive::File &f, const std::vector<double> &profile,
    const phmap::btree_map<hictk::Chromosome, double> &scaling_factors) {
  auto grp = f.createGroup("profile");
  if (!profile.empty()) {
    HighFive::DataSetCreateProps props{};
    props.add(HighFive::Chunking({profile.size()}));
    props.add(HighFive::Deflate(9));
    grp.createDataSet("values", profile, props);

    std::vector<double> scaling_factors_flat(scaling_factors.size());
    std::transform(scaling_factors.begin(), scaling_factors.end(), scaling_factors_flat.begin(),
                   [](const auto &kv) { return kv.second; });
    grp.createDataSet("scaling-factors", scaling_factors_flat);
  }
}

template <typename File>
inline void ExpectedValues<File>::serialize_trans_profiles(
    HighFive::File &f, const phmap::btree_map<EVTKey, double> &nnz_avg_values) {
  auto grp = f.createGroup("avg-values");
  if (nnz_avg_values.empty()) {
    return;
  }

  std::vector<std::string> chrom1{};
  std::vector<std::string> chrom2{};
  std::vector<double> values{};

  for (const auto &[cp, value] : nnz_avg_values) {
    chrom1.emplace_back(std::string{cp.first.name()});   // NOLINT
    chrom2.emplace_back(std::string{cp.second.name()});  // NOLINT
    values.emplace_back(value);
  }

  grp.createDataSet("chrom1", chrom1);
  grp.createDataSet("chrom2", chrom2);
  grp.createDataSet("value", values);
}

template <typename File>
inline hictk::Reference ExpectedValues<File>::deserialize_chromosomes(const HighFive::File &f) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  f.getDataSet("chroms/name").read(chrom_names);
  f.getDataSet("chroms/length").read(chrom_sizes);

  assert(chrom_names.size() == chrom_sizes.size());
  return {chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};
}

template <typename File>
inline std::pair<const std::vector<double>, const phmap::btree_map<hictk::Chromosome, double>>
ExpectedValues<File>::deserialize_cis_profiles(const HighFive::File &f) {
  std::vector<double> weights{};
  if (f.exist("profile/values")) {
    f.getDataSet("profile/values").read(weights);
  }

  phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
  if (f.exist("profile/scaling-factors")) {
    const auto chroms = deserialize_chromosomes(f);
    std::vector<double> values{};
    f.getDataSet("profile/scaling-factors").read(values);

    assert(chroms.size() == values.size());
    for (std::uint32_t i = 0; i < chroms.size(); ++i) {
      scaling_factors.emplace(chroms[i], values[i]);
    }
  }

  return std::make_pair(weights, scaling_factors);
}

template <typename File>
inline auto ExpectedValues<File>::deserialize_trans_profiles(const HighFive::File &f)
    -> phmap::btree_map<EVTKey, double> {
  if (!f.exist("avg-values/value")) {
    return {};
  }

  const auto chroms = deserialize_chromosomes(f);
  std::vector<std::string> chrom1{};
  std::vector<std::string> chrom2{};
  std::vector<double> values{};

  f.getDataSet("avg-values/chrom1").read(chrom1);
  f.getDataSet("avg-values/chrom2").read(chrom2);
  f.getDataSet("avg-values/value").read(values);

  assert(chrom1.size() == chrom2.size());
  assert(chrom1.size() == values.size());

  phmap::btree_map<EVTKey, double> expected_values{};
  for (std::size_t i = 0; i < values.size(); ++i) {
    expected_values.emplace(EVTKey{chroms.at(chrom1[i]), chroms.at(chrom2[i])}, values[i]);
  }
  return expected_values;
}

}  // namespace nchg
