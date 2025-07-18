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

#include "nchg/expected_values.hpp"

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
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

#include "nchg/chromosome_pairs_generator.hpp"
#include "nchg/concepts.hpp"
#include "nchg/expected_matrix.hpp"
#include "nchg/expected_values_aggregator.hpp"
#include "nchg/mad_max_filter.hpp"

namespace nchg {

ExpectedValues::ExpectedValues(
    std::shared_ptr<const hictk::File> file, const Params &params_,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask)
    : _fp(std::move(file)),
      _resolution(!!_fp ? _fp->resolution() : std::uint32_t{0}),
      _mad_max(params_.mad_max),
      _min_delta(params_.min_delta),
      _max_delta(params_.max_delta),
      _bin_aggregation_possible_distances_cutoff(params_.bin_aggregation_possible_distances_cutoff),
      _bin_aggregation_observed_distances_cutoff(params_.bin_aggregation_observed_distances_cutoff),
      _interpolate(params_.interpolate),
      _interpolation_qtile(params_.interpolation_qtile),
      _interpolation_window_size(params_.interpolation_window_size),
      _seeded(!bin_mask.empty()) {
  if (_mad_max < 0 || !std::isfinite(_mad_max)) {
    throw std::logic_error("mad_max should be a non-negative value");
  }
  if (_min_delta > _max_delta) {
    throw std::logic_error("min_delta should be less than or equal to max_delta");
  }

  if (_interpolation_qtile < 0 || _interpolation_qtile > 1) {
    throw std::logic_error("interpolation_qtile should be between 0 and 1");
  }

  if (_fp) {
    compute_expected_values_cis(bin_mask);
    compute_expected_values_trans(bin_mask);
  }
}

ExpectedValues ExpectedValues::cis_only(
    std::shared_ptr<const hictk::File> file, const Params &params_,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask) {
  ExpectedValues ev(nullptr, params_);
  ev._fp = std::move(file);
  ev._resolution = ev._fp->resolution();
  ev._seeded = !bin_mask.empty();
  if (ev._fp) {
    ev.compute_expected_values_cis(bin_mask);
  }
  return ev;
}

ExpectedValues ExpectedValues::trans_only(
    std::shared_ptr<const hictk::File> file, const Params &params_,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask) {
  ExpectedValues ev(nullptr, {.mad_max = params_.mad_max,
                              .min_delta = 0,
                              .max_delta = std::numeric_limits<std::uint64_t>::max(),
                              .bin_aggregation_possible_distances_cutoff = 0,
                              .bin_aggregation_observed_distances_cutoff = 0,
                              .interpolate = false,
                              .interpolation_qtile = 0,
                              .interpolation_window_size = 0});
  ev._fp = std::move(file);
  ev._resolution = ev._fp->resolution();
  ev._seeded = !bin_mask.empty();
  if (ev._fp) {
    ev.compute_expected_values_trans(bin_mask);
  }
  return ev;
}

ExpectedValues ExpectedValues::chromosome_pair(
    std::shared_ptr<const hictk::File> file, const hictk::Chromosome &chrom1,
    const hictk::Chromosome &chrom2, const Params &params,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask) {
  SPDLOG_INFO("[{}:{}]: computing expected values...", chrom1.name(), chrom2.name());

  ExpectedValues ev(nullptr, params);
  ev._fp = std::move(file);
  ev._resolution = ev._fp->resolution();
  ev._seeded = !bin_mask.empty();
  if (chrom1 == chrom2) {
    ev.compute_expected_values_cis(bin_mask);
    return ev;
  }

  std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom1.name(), chrom2.name());
        const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<N>(),
                                                          sel.template end<N>(), fp.bins_ptr());

        if (!ev._bin_masks.contains(std::make_pair(chrom1, chrom2))) {
          ev.add_bin_mask(chrom1, chrom2,
                          mad_max_filtering(std::ranges::subrange(jsel.begin(), jsel.end()), chrom1,
                                            chrom2, ev._resolution, ev._mad_max),
                          bin_mask);
        }

        const ExpectedMatrixStats em(jsel, chrom1, chrom2, ev._fp->bins(), std::vector<N>{}, 0,
                                     *ev.bin_mask(chrom1, chrom2).first,
                                     *ev.bin_mask(chrom1, chrom2).second,
                                     std::numeric_limits<std::uint64_t>::max());

        ev._expected_values_trans.emplace(std::make_pair(chrom1, chrom2), em.nnz_avg());
      },
      ev._fp->get());

  return ev;
}

ExpectedValues ExpectedValues::deserialize(const std::filesystem::path &path) {
  ExpectedValues ev{nullptr};
  HighFive::File f(path.string());

  ev._resolution = f.getAttribute("resolution").read<std::uint32_t>();

  try {
    ev._seeded = f.getGroup("bin-masks").getAttribute("seeded").read<bool>();
  } catch (const HighFive::Exception &) {
    ev._seeded = false;
  }

  const auto params = deserialize_attributes(f);
  ev._mad_max = params.mad_max;
  ev._min_delta = params.min_delta;
  ev._max_delta = params.max_delta;
  ev._bin_aggregation_possible_distances_cutoff = params.bin_aggregation_possible_distances_cutoff;
  ev._bin_aggregation_observed_distances_cutoff = params.bin_aggregation_observed_distances_cutoff;
  ev._interpolate = params.interpolate;
  ev._interpolation_qtile = params.interpolation_qtile;
  ev._interpolation_window_size = params.interpolation_window_size;

  ev._bin_masks = deserialize_bin_masks(f);
  auto [weights, scaling_factors] = deserialize_cis_profiles(f);
  ev._expected_weights = std::move(weights);
  ev._expected_scaling_factors = std::move(scaling_factors);
  ev._expected_values_trans = deserialize_trans_profiles(f);

  return ev;
}

std::uint32_t ExpectedValues::resolution() const noexcept { return _resolution; }

auto ExpectedValues::params() const noexcept -> Params {
  return {.mad_max = _mad_max,
          .min_delta = _min_delta,
          .max_delta = _max_delta,
          .bin_aggregation_possible_distances_cutoff = _bin_aggregation_possible_distances_cutoff,
          .bin_aggregation_observed_distances_cutoff = _bin_aggregation_observed_distances_cutoff,
          .interpolate = _interpolate,
          .interpolation_qtile = _interpolation_qtile,
          .interpolation_window_size = _interpolation_window_size};
}

std::shared_ptr<const std::vector<bool>> ExpectedValues::bin_mask(
    const hictk::Chromosome &chrom) const {
  return _bin_masks.at(std::make_pair(chrom, chrom)).first;
}

std::pair<std::shared_ptr<const std::vector<bool>>, std::shared_ptr<const std::vector<bool>>>
ExpectedValues::bin_mask(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const {
  return _bin_masks.at(std::make_pair(chrom1, chrom2));
}

std::vector<double> ExpectedValues::expected_values(const hictk::Chromosome &chrom,
                                                    bool rescale) const {
  if (_expected_weights.empty()) {
    throw std::out_of_range(
        fmt::format("expected values for \"{}\" are not available: out of range", chrom.name()));
  }

  const auto num_bins = (chrom.size() + _resolution - 1) / _resolution;
  assert(num_bins <= _expected_weights.size());

  std::vector weights(_expected_weights.begin(),
                      _expected_weights.begin() + static_cast<std::ptrdiff_t>(num_bins));

  const auto sf = rescale ? _expected_scaling_factors.at(chrom) : 1.0;
  if (sf == 1) {
    return weights;
  }

  std::ranges::transform(weights, weights.begin(), [&](auto n) {
    n /= sf;
    if (std::isfinite(n)) {
      return n;
    }
    return 0.0;
  });
  return weights;
}

double ExpectedValues::expected_value(const hictk::Chromosome &chrom1,
                                      const hictk::Chromosome &chrom2) const {
  return _expected_values_trans.at(ChromPair{chrom1, chrom2});
}

const phmap::btree_map<hictk::Chromosome, double> &ExpectedValues::scaling_factors()
    const noexcept {
  return _expected_scaling_factors;
}

double ExpectedValues::scaling_factor(const hictk::Chromosome &chrom) const {
  return scaling_factors().at(chrom);
}

ExpectedMatrixStats ExpectedValues::expected_matrix(const hictk::Chromosome &chrom) const {
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  return std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom.name());
        const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<N>(),
                                                          sel.template end<N>(), fp.bins_ptr());
        return expected_matrix(chrom, _fp->bins(), jsel);
      },
      _fp->get());
}

ExpectedMatrixStats ExpectedValues::expected_matrix(const hictk::Chromosome &chrom1,
                                                    const hictk::Chromosome &chrom2) const {
  if (chrom1 == chrom2) {
    return expected_matrix(chrom1);
  }

  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  return std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom1.name(), chrom2.name());
        const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<N>(),
                                                          sel.template end<N>(), fp.bins_ptr());
        return expected_matrix(chrom1, chrom2, fp.bins(), jsel);
      },
      _fp->get());
}

void ExpectedValues::serialize(const std::filesystem::path &path) const {
  SPDLOG_INFO("writing expected value profiles to {}...", path);
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  HighFive::File f(path.string(), HighFive::File::Create);

  const auto source_file = std::filesystem::path{_fp->path()}.filename().string();
  f.createAttribute("source-file", source_file);
  f.createAttribute("resolution", _resolution);
  serialize_attributes(f, params());

  serialize_chromosomes(f, _fp->chromosomes());
  serialize_bin_masks(f, _bin_masks, _seeded);
  serialize_cis_profiles(f, _expected_weights, _expected_scaling_factors);
  serialize_trans_profiles(f, _expected_values_trans);
}

bool ExpectedValues::cis_only() const noexcept {
  return !_expected_weights.empty() && _expected_values_trans.empty();
}

bool ExpectedValues::trans_only() const noexcept {
  return _expected_weights.empty() && !_expected_values_trans.empty();
}

void ExpectedValues::compute_expected_values_cis(
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask_seed) {
  SPDLOG_INFO("initializing expected matrix weights from cis interactions...");
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  ExpectedValuesAggregator aggr(_fp->bins_ptr());

  const auto no_pixels_processed = std::visit(
      [&](const auto &f) {
        init_bin_masks(f, bin_mask_seed);

        const auto [selectors, merger] = init_pixel_merger_cis(f);
        if (selectors.empty()) {
          return true;
        }

        const hictk::transformers::JoinGenomicCoords mjsel(merger.begin(), merger.end(),
                                                           f.bins_ptr());

        for (const hictk::Pixel<N> &p : mjsel) {
          const auto &mask = *bin_mask(p.coords.bin1.chrom());
          const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
          const auto bin1_id = p.coords.bin1.rel_id();
          const auto bin2_id = p.coords.bin2.rel_id();
          const auto bin1_masked = !mask.empty() && mask[bin1_id];
          const auto bin2_masked = !mask.empty() && mask[bin2_id];
          if (delta >= _min_delta && delta < _max_delta && !bin1_masked && !bin2_masked)
              [[likely]] {
            aggr.add(p);
          }
        }

        return false;
      },
      _fp->get());

  if (no_pixels_processed) {
    _expected_weights.resize(_fp->bins().size(), 0);
    for (const auto &chrom : _fp->chromosomes()) {
      _expected_scaling_factors.emplace(chrom, 1.0);
    }
    return;
  }

  aggr.compute_density(_bin_aggregation_possible_distances_cutoff,
                       _bin_aggregation_observed_distances_cutoff, _interpolate,
                       _interpolation_qtile, _interpolation_window_size);

  _expected_scaling_factors = aggr.scaling_factors();
  for (const auto &chrom : _fp->chromosomes()) {
    if (!_expected_scaling_factors.contains(chrom)) {
      _expected_scaling_factors.emplace(chrom, std::numeric_limits<double>::quiet_NaN());
    }
    SPDLOG_DEBUG("[{}]: scaling_factor={}", chrom.name(), _expected_scaling_factors.at(chrom));
  }

  _expected_weights = aggr.weights();
  const auto &chrom = _fp->chromosomes().longest_chromosome();
  const auto num_bins = (chrom.size() + _resolution - 1) / _resolution;
  _expected_weights.resize(num_bins, std::numeric_limits<double>::quiet_NaN());
  SPDLOG_INFO("finished computing the expected value profile for cis interactions");
}

void ExpectedValues::compute_expected_values_trans(
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask_seed) {
  if (!_fp) {
    throw std::logic_error("ExpectedValues::expected_matrix() was called on a null file");
  }

  std::visit(
      [&](const auto &f) {
        for (const auto &[chrom1, chrom2] :
             generate_chromosome_pairs_upper_triangle(f.chromosomes(), false)) {
          if (chrom1 == chrom2) {
            continue;
          }

          const auto sel = f.fetch(chrom1.name(), chrom2.name());
          const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<N>(),
                                                            sel.template end<N>(), f.bins_ptr());

          if (!_bin_masks.contains(std::make_pair(chrom1, chrom2))) {
            if (sel.empty()) {
              SPDLOG_DEBUG(
                  "[{}:{}]: no interactions found: masking the entire chromosome-chromosome "
                  "matrix!",
                  chrom1.name(), chrom2.name());
              const auto num_bins1 = (chrom1.size() + f.resolution() - 1) / f.resolution();
              const auto num_bins2 = (chrom2.size() + f.resolution() - 1) / f.resolution();
              auto masks =
                  std::make_pair(std::vector(num_bins1, true), std::vector(num_bins2, true));
              add_bin_mask(chrom1, chrom2, std::move(masks), bin_mask_seed);
            } else {
              SPDLOG_INFO("[{}:{}]: begin computing expected value...", chrom1.name(),
                          chrom2.name());
              add_bin_mask(chrom1, chrom2,
                           mad_max_filtering(jsel, chrom1, chrom2, _resolution, _mad_max),
                           bin_mask_seed);
            }
          }

          const ExpectedMatrixStats em(
              jsel, chrom1, chrom2, f.bins(), std::vector<N>{}, 0, *bin_mask(chrom1, chrom2).first,
              *bin_mask(chrom1, chrom2).second, std::numeric_limits<std::uint64_t>::max());
          SPDLOG_DEBUG("[{}:{}]: expected_value={}", chrom1.name(), chrom2.name(), em.nnz_avg());
          _expected_values_trans.emplace(std::make_pair(chrom1, chrom2), em.nnz_avg());
        }
      },
      _fp->get());
}

void ExpectedValues::add_bin_mask(
    const hictk::Chromosome &chrom, std::vector<bool> mask,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask_seed) {
  const auto key = std::make_pair(chrom, chrom);
  if (_bin_masks.contains(key)) {
    return;
  }
  const auto value = std::make_shared<std::vector<bool>>(std::move(mask));
  const auto match = bin_mask_seed.find(chrom);
  if (match != bin_mask_seed.end()) {
    merge_bin_masks(*value, match->second);
  }
  _bin_masks.emplace(key, std::make_pair(std::const_pointer_cast<const std::vector<bool>>(value),
                                         std::const_pointer_cast<const std::vector<bool>>(value)));
}

void ExpectedValues::add_bin_mask(
    const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
    std::pair<std::vector<bool>, std::vector<bool>> masks,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask_seed) {
  const auto key = std::make_pair(chrom1, chrom2);
  if (_bin_masks.contains(key)) {
    return;
  }
  auto value1 = std::make_shared<std::vector<bool>>(std::move(masks.first));
  auto value2 = std::make_shared<std::vector<bool>>(std::move(masks.second));

  auto match1 = bin_mask_seed.find(chrom1);
  if (match1 != bin_mask_seed.end()) {
    merge_bin_masks(*value1, match1->second);
  }
  auto match2 = bin_mask_seed.find(chrom2);
  if (match2 != bin_mask_seed.end()) {
    merge_bin_masks(*value2, match2->second);
  }

  _bin_masks.emplace(
      key, std::make_pair(std::const_pointer_cast<const std::vector<bool>>(std::move(value1)),
                          std::const_pointer_cast<const std::vector<bool>>(std::move(value2))));
}

void ExpectedValues::serialize_attributes(HighFive::File &f, const Params &params) {
  assert(params.min_delta <= params.max_delta);
  f.createAttribute("mad_max", params.mad_max);
  f.createAttribute("min_delta", params.min_delta);
  f.createAttribute("max_delta", params.max_delta);
  f.createAttribute("bin_aggregation_possible_distances_cutoff",
                    params.bin_aggregation_possible_distances_cutoff);
  f.createAttribute("bin_aggregation_observed_distances_cutoff",
                    params.bin_aggregation_observed_distances_cutoff);
  f.createAttribute("interpolate", params.interpolate);
  f.createAttribute("interpolation_qtile", params.interpolation_qtile);
  f.createAttribute("interpolation_window_size", params.interpolation_window_size);
}

void ExpectedValues::serialize_chromosomes(HighFive::File &f, const hictk::Reference &chroms) {
  auto grp = f.createGroup("chroms");

  std::vector<std::string> chromosome_names(chroms.size());
  std::ranges::transform(chroms, chromosome_names.begin(),
                         [](const hictk::Chromosome &c) { return std::string{c.name()}; });
  grp.createDataSet("name", chromosome_names);

  std::vector<std::uint32_t> chromosome_sizes(chroms.size());
  std::ranges::transform(chroms, chromosome_sizes.begin(),
                         [](const hictk::Chromosome &c) { return c.size(); });
  grp.createDataSet("length", chromosome_sizes);
}

void ExpectedValues::serialize_bin_masks(
    HighFive::File &f, const phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>> &bin_masks,
    bool seeded) {
  auto grp = f.createGroup("bin-masks");
  grp.createAttribute("seeded", seeded);
  if (bin_masks.empty()) {
    return;
  }

  std::vector<std::string> chrom1{};
  std::vector<std::string> chrom2{};
  std::vector<std::size_t> offsets1{};
  std::vector<std::size_t> offsets2{};
  std::vector<bool> values1{};
  std::vector<bool> values2{};

  for (const auto &[k, v] : bin_masks) {
    chrom1.emplace_back(k.first.name());
    chrom2.emplace_back(k.second.name());

    offsets1.emplace_back(values1.size());
    values1.insert(values1.end(), v.first->begin(), v.first->end());

    offsets2.emplace_back(values2.size());
    values2.insert(values2.end(), v.second->begin(), v.second->end());
  }

  offsets1.emplace_back(values1.size());
  offsets2.emplace_back(values2.size());

  grp.createDataSet("chrom1", chrom1);
  grp.createDataSet("chrom2", chrom2);
  grp.createDataSet("offsets1", offsets1);
  grp.createDataSet("offsets2", offsets2);

  HighFive::DataSetCreateProps props1{};
  props1.add(HighFive::Chunking({values1.size()}));
  props1.add(HighFive::Deflate(9));

  grp.createDataSet("values1", values1, props1);

  HighFive::DataSetCreateProps props2{};
  props2.add(HighFive::Chunking({values2.size()}));
  props2.add(HighFive::Deflate(9));
  grp.createDataSet("values2", values2, props2);
}

void ExpectedValues::serialize_cis_profiles(
    HighFive::File &f, const std::vector<double> &profile,
    const phmap::btree_map<hictk::Chromosome, double> &scaling_factors) {
  auto grp = f.createGroup("profile");
  if (!profile.empty()) {
    HighFive::DataSetCreateProps props{};
    props.add(HighFive::Chunking({profile.size()}));
    props.add(HighFive::Deflate(9));
    grp.createDataSet("values", profile, props);

    std::vector<double> scaling_factors_flat(scaling_factors.size());
    std::ranges::transform(scaling_factors, scaling_factors_flat.begin(),
                           [](const auto &kv) { return kv.second; });
    grp.createDataSet("scaling-factors", scaling_factors_flat);
  }
}

void ExpectedValues::serialize_trans_profiles(
    HighFive::File &f, const phmap::btree_map<ChromPair, double> &nnz_avg_values) {
  auto grp = f.createGroup("avg-values");
  if (nnz_avg_values.empty()) {
    return;
  }

  std::vector<std::string> chrom1{};
  std::vector<std::string> chrom2{};
  std::vector<double> values{};

  for (const auto &[cp, value] : nnz_avg_values) {
    chrom1.emplace_back(cp.first.name());
    chrom2.emplace_back(cp.second.name());
    values.emplace_back(value);
  }

  grp.createDataSet("chrom1", chrom1);
  grp.createDataSet("chrom2", chrom2);
  grp.createDataSet("value", values);
}

auto ExpectedValues::deserialize_attributes(const HighFive::File &f) -> Params {
  return {.mad_max = f.getAttribute("mad_max").read<double>(),
          .min_delta = f.getAttribute("min_delta").read<std::uint64_t>(),
          .max_delta = f.getAttribute("max_delta").read<std::uint64_t>(),
          .bin_aggregation_possible_distances_cutoff =
              f.getAttribute("bin_aggregation_possible_distances_cutoff").read<double>(),
          .bin_aggregation_observed_distances_cutoff =
              f.getAttribute("bin_aggregation_observed_distances_cutoff").read<double>(),
          .interpolate = f.getAttribute("interpolate").read<bool>(),
          .interpolation_qtile = f.getAttribute("interpolation_qtile").read<double>(),
          .interpolation_window_size =
              f.getAttribute("interpolation_window_size").read<std::uint32_t>()};
}

hictk::Reference ExpectedValues::deserialize_chromosomes(const HighFive::File &f) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  f.getDataSet("chroms/name").read(chrom_names);
  f.getDataSet("chroms/length").read(chrom_sizes);

  assert(chrom_names.size() == chrom_sizes.size());
  return {chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};
}

auto ExpectedValues::deserialize_bin_masks(HighFive::File &f)
    -> phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>> {
  const auto grp = f.getGroup("bin-masks");
  if (!grp.exist("chrom1")) {
    return {};
  }

  const auto chroms = deserialize_chromosomes(f);

  const auto chrom1_names = grp.getDataSet("chrom1").read<std::vector<std::string>>();
  const auto chrom2_names = grp.getDataSet("chrom2").read<std::vector<std::string>>();
  const auto offsets1 = grp.getDataSet("offsets1").read<std::vector<std::size_t>>();
  const auto offsets2 = grp.getDataSet("offsets2").read<std::vector<std::size_t>>();
  const auto values1 = grp.getDataSet("values1").read<std::vector<bool>>();
  const auto values2 = grp.getDataSet("values2").read<std::vector<bool>>();

  phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>> buffer{};
  for (std::size_t i = 0; i < chrom1_names.size(); ++i) {
    const auto &chrom1 = chroms.at(chrom1_names.at(i));
    const auto &chrom2 = chroms.at(chrom2_names.at(i));

    const auto j0 = offsets1.at(i);
    const auto j1 = offsets1.at(i + 1);

    const auto k0 = offsets2.at(i);
    const auto k1 = offsets2.at(i + 1);

    std::vector<bool> mask1(j1 - j0);
    std::vector<bool> mask2{};

    std::copy(values1.begin() + static_cast<std::ptrdiff_t>(j0),
              values1.begin() + static_cast<std::ptrdiff_t>(j1), mask1.begin());
    if (chrom1 != chrom2) {
      mask2.resize(k1 - k0);
      std::copy(values2.begin() + static_cast<std::ptrdiff_t>(k0),
                values2.begin() + static_cast<std::ptrdiff_t>(k1), mask2.begin());
    }

    auto mask1_ptr = std::make_shared<const std::vector<bool>>(std::move(mask1));
    auto mask2_ptr = mask2.empty() ? mask1_ptr : std::make_shared<const std::vector<bool>>(mask2);
    buffer.emplace(std::make_pair(chrom1, chrom2),
                   std::make_pair(std::move(mask1_ptr), std::move(mask2_ptr)));
  }

  return buffer;
}

std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
ExpectedValues::deserialize_cis_profiles(const HighFive::File &f) {
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

  return {weights, scaling_factors};
}

auto ExpectedValues::deserialize_trans_profiles(const HighFive::File &f)
    -> phmap::btree_map<ChromPair, double> {
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

  phmap::btree_map<ChromPair, double> expected_values{};
  for (std::size_t i = 0; i < values.size(); ++i) {
    expected_values.emplace(ChromPair{chroms.at(chrom1[i]), chroms.at(chrom2[i])}, values[i]);
  }
  return expected_values;
}

void ExpectedValues::merge_bin_masks(std::vector<bool> &mask1, const std::vector<bool> &mask2) {
  if (mask1.size() != mask2.size()) {
    throw std::runtime_error(fmt::format("bin mask shape mismatch: expected shape {}, found {}",
                                         mask1.size(), mask2.size()));
  }
  for (std::size_t i = 0; i < mask1.size(); ++i) {
    mask1[i] = mask1[i] || mask2[i];
  }
}

}  // namespace nchg
