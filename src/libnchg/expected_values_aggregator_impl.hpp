// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>

#pragma once

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"
#include "nchg/common.hpp"

namespace nchg {

inline ExpectedValuesAggregator::ExpectedValuesAggregator(
    std::shared_ptr<const hictk::BinTable> bins)
    : _bins(std::move(bins)) {
  SPDLOG_INFO(FMT_STRING("[{} bp] initializing expected value vector"), _bins->resolution());
  std::uint32_t max_length = 0;
  for (const auto &chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }

    max_length = std::max(max_length, chrom.size());

    _num_bins_gw += chrom.size();
  }

  const auto bin_size = _bins->resolution();
  // round down to mimick HiCTools' behavior
  const auto max_n_bins = max_length / bin_size;
  _possible_distances.resize(max_n_bins, 0.0);
  _observed_distances.resize(max_n_bins, 0.0);
  _weights.resize(max_n_bins, 0.0);
}

template <typename N>
inline void ExpectedValuesAggregator::add(const hictk::ThinPixel<N> &p) {
  add(hictk::Pixel<N>{*_bins, p});
}

template <typename N>
inline void ExpectedValuesAggregator::add(const hictk::Pixel<N> &p) {
  const auto count = conditional_static_cast<double>(p.count);
  if (std::isnan(count)) {
    return;
  }

  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();

  if (p.coords.is_intra()) {
    at(chrom1) += count;
    const auto i = p.coords.bin2.id() - p.coords.bin1.id();
    // skip last bin in chromosome if chromosome size is not a multiple of bin size
    // this is done to mimick HiCTools' behavior
    if (i < _observed_distances.size()) {
      _observed_distances[i] += count;
    }
  } else {
    at(chrom1, chrom2) += count;
  }
}

inline void ExpectedValuesAggregator::compute_density(
    double bin_aggregation_possible_distances_cutoff,
    double bin_aggregation_observed_distances_cutoff, bool intepolate_weights,
    double interpolation_qtile, std::uint32_t interpolation_window_size) {
  SPDLOG_INFO(FMT_STRING("[{} bp] computing expected vector density"), _bins->resolution());
  init_possible_distances();
  compute_density_cis(bin_aggregation_possible_distances_cutoff,
                      bin_aggregation_observed_distances_cutoff, intepolate_weights,
                      interpolation_qtile, interpolation_window_size);
  compute_density_trans();
}

inline const std::vector<double> &ExpectedValuesAggregator::weights() const noexcept {
  return _weights;
}

inline std::vector<double> ExpectedValuesAggregator::weights(const hictk::Chromosome &chrom,
                                                             bool rescale) const {
  const auto num_bins = chrom.size() / _bins->resolution();

  std::vector<double> w{_weights.begin(), _weights.begin() + static_cast<std::ptrdiff_t>(num_bins)};
  if (!rescale) {
    return w;
  }

  const auto sf = scaling_factor(chrom);
  std::transform(w.begin(), w.end(), w.begin(), [&](const auto n) { return n / sf; });

  return w;
}

inline double ExpectedValuesAggregator::scaling_factor(const hictk::Chromosome &chrom) const {
  return _scaling_factors.at(chrom);
}

inline const phmap::btree_map<hictk::Chromosome, double> &
ExpectedValuesAggregator::scaling_factors() const noexcept {
  return _scaling_factors;
}

inline void ExpectedValuesAggregator::init_possible_distances() {
  const auto bin_size = _bins->resolution();

  for (const auto &[chrom, _] : _cis_sum) {
    if (chrom.is_all()) {
      continue;
    }
    const auto n_bins = chrom.size() / bin_size;
    for (std::uint32_t i = 0; i < n_bins; ++i) {
      _possible_distances[i] += n_bins - i;
    }
  }
}

inline void ExpectedValuesAggregator::compute_density_cis(
    double bin_aggregation_possible_distances_cutoff,
    double bin_aggregation_observed_distances_cutoff, bool interpolate_weights,
    double interpolation_qtile, std::uint32_t interpolation_window_size) {
  if (_observed_distances.empty()) {
    return;
  }
  const auto last_valid_bin_idx = aggregate_low_cov_bins(bin_aggregation_possible_distances_cutoff,
                                                         bin_aggregation_observed_distances_cutoff);

  assert(_weights.size() == _observed_distances.size());
  for (std::size_t i = 0; i < _weights.size(); i++) {
    _weights[i] = _observed_distances[i] == 0 ? std::numeric_limits<double>::quiet_NaN()
                                              : _observed_distances[i] / _possible_distances[i];
  }

  if (last_valid_bin_idx != _possible_distances.size()) {
    const auto it1 = _weights.begin() + static_cast<std::ptrdiff_t>(last_valid_bin_idx);
    const auto it2 = _weights.end();
    std::fill(it1, it2, std::numeric_limits<double>::quiet_NaN());
  }

  if (interpolate_weights) {
    correct_outliers(interpolation_qtile, interpolation_window_size);
  }

  _scaling_factors = compute_scaling_factors();
}

inline std::size_t ExpectedValuesAggregator::aggregate_low_cov_bins(
    double possible_distances_cutoff, double observed_distances_cutoff) noexcept {
  std::size_t i = 0;
  while (_observed_distances[i] == 0 && i < _weights.size()) {
    ++i;
  }

  auto last_nnz = std::find_if(_observed_distances.rbegin(), _observed_distances.rend(),
                               [](const auto n) { return n != 0; });
  if (last_nnz != _observed_distances.rbegin()) {
    --last_nnz;
  }

  std::size_t i1{};
  for (; i < _observed_distances.size(); ++i) {
    auto pos_sum = _possible_distances[i];
    auto obs_sum = _observed_distances[i];
    std::size_t j = 1;
    std::size_t n = 1;
    while ((pos_sum < possible_distances_cutoff || obs_sum < observed_distances_cutoff) && i >= j &&
           i + j < _weights.size()) {
      pos_sum += _possible_distances[i - j] + _possible_distances[i + j];
      obs_sum += _observed_distances[i - j] + _observed_distances[i + j];
      n += 2;
      ++j;
    }

    if (pos_sum >= possible_distances_cutoff || obs_sum >= observed_distances_cutoff) {
      i1 = i + j;
      _possible_distances[i] = pos_sum / static_cast<double>(n);
      _observed_distances[i] = obs_sum / static_cast<double>(n);
    }
  }

  i1 = std::min(i1, static_cast<std::size_t>(std::distance(last_nnz, _observed_distances.rend())));

  auto it1 = _observed_distances.begin() + static_cast<std::ptrdiff_t>(i1);
  auto it2 = _observed_distances.end();
  std::fill(it1, it2, std::numeric_limits<double>::quiet_NaN());

  return i1;
}

inline phmap::btree_map<hictk::Chromosome, double>
ExpectedValuesAggregator::compute_scaling_factors() const {
  phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
  for (const auto &[chrom, _] : _cis_sum) {
    if (chrom.is_all()) {
      continue;
    }

    const auto num_chrom_bins =
        std::min(_weights.size(), chrom.size() / static_cast<std::size_t>(_bins->resolution()));
    auto expected_count = 0.0;
    for (std::size_t i = 0; i < num_chrom_bins; ++i) {
      const auto v = _weights[i];
      if (std::isfinite(v)) {
        expected_count += (static_cast<double>(num_chrom_bins) - static_cast<double>(i)) * v;
      }
    }

    const auto observed_count = _cis_sum.at(chrom);
    const auto f = expected_count / observed_count;
    scaling_factors.emplace(chrom, f);
  }
  return scaling_factors;
}

namespace internal {

[[nodiscard]] inline double median(std::vector<double> &v) {
  if (v.empty()) {
    throw std::runtime_error("median was called on an empty vector");
  }

  const auto size = static_cast<std::ptrdiff_t>(v.size());
  auto first = v.begin();
  auto mid = first + (size / 2);
  auto last = v.end();

  std::nth_element(first, mid, last);

  if (size % 2 != 0) {
    return *mid;
  }

  const auto n1 = *mid;
  std::nth_element(first, --mid, last);
  const auto n2 = *mid;

  return (n1 + n2) / 2;
}

[[nodiscard]] inline double approx_quantile(std::vector<double> &v, double quantile) {
  auto nth = v.begin() + static_cast<std::ptrdiff_t>(quantile * static_cast<double>(v.size() - 1));
  std::nth_element(v.begin(), nth, v.end());
  return *nth;
};

}  // namespace internal

inline void ExpectedValuesAggregator::correct_outliers(double quantile, std::uint32_t window_size) {
  assert(quantile >= 0);
  assert(quantile <= 1);

  const auto window_size_bins = (window_size + _bins->resolution() - 1) / _bins->resolution();
  if (window_size_bins < 2) {
    return;
  }

  const auto window_arm = window_size_bins / 2;

  std::vector<double> buffer{};
  for (std::size_t i = window_arm; i + window_arm < _weights.size(); ++i) {
    if (!std::isfinite(_weights[i])) {
      continue;
    }
    const auto i0 = i - window_arm;
    const auto i1 = i + window_arm + 1;
    buffer.clear();
    std::copy(_weights.begin() + static_cast<std::ptrdiff_t>(i0),
              _weights.begin() + static_cast<std::ptrdiff_t>(i1), std::back_inserter(buffer));
    const auto thresh = internal::approx_quantile(buffer, quantile);
    if (_weights[i] > thresh) {
      const auto n = internal::median(buffer);
      if (std::isfinite(n)) {
        _weights[i] = n;
      }
    }
  }
}

inline void ExpectedValuesAggregator::compute_density_trans() {
  for (auto &[k, v] : _trans_sum) {
    // We round-down to match HiCTools behavior
    const auto num_bins1 = k.first.size() / _bins->resolution();
    const auto num_bins2 = k.second.size() / _bins->resolution();
    const auto num_pixels = num_bins1 * num_bins2;
    v = num_pixels != 0 ? v / static_cast<double>(num_pixels) : 0.0;
  }
}

inline double ExpectedValuesAggregator::at(const hictk::Chromosome &chrom) const {
  return _cis_sum.at(chrom);
}

inline double ExpectedValuesAggregator::at(const hictk::Chromosome &chrom1,
                                           const hictk::Chromosome &chrom2) const {
  return _trans_sum.at(std::make_pair(chrom1, chrom2));
}

inline double &ExpectedValuesAggregator::at(const hictk::Chromosome &chrom) {
  auto [it, _] = _cis_sum.try_emplace(chrom, 0.0);
  return it->second;
}

inline double &ExpectedValuesAggregator::at(const hictk::Chromosome &chrom1,
                                            const hictk::Chromosome &chrom2) {
  auto [it, _] = _trans_sum.try_emplace(std::make_pair(chrom1, chrom2), 0.0);
  return it->second;
}

inline const hictk::Reference &ExpectedValuesAggregator::chromosomes() const noexcept {
  assert(_bins);
  return _bins->chromosomes();
}

}  // namespace nchg
