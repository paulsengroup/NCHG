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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/expected_values_aggregator.hpp>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/concepts.hpp"
#include "nchg/matrix_stats.hpp"

namespace nchg {

namespace internal {
template <typename N>
constexpr bool is_close(N n1, N n2, double rtol = 1.0e-6, double atol = 0) noexcept {
  assert(rtol >= 0 && rtol <= 1);
  if constexpr (std::is_integral_v<N>) {
    return n1 == n2;
  } else {
    if (n1 == n2) {
      return true;
    }
    if (std::isnan(n1)) {
      return std::isnan(n2);
    }

    // https://peps.python.org/pep-0485/
    const auto diff = std::abs(n1 - n2);
    return (diff <= abs(rtol * n2)) || (diff <= atol);
  }
}
}  // namespace internal

template <typename Pixels>
  requires PixelRange<Pixels>
inline ExpectedMatrixStats::ExpectedMatrixStats(const Pixels &pixels, hictk::Chromosome chrom1,
                                                hictk::Chromosome chrom2, hictk::BinTable bins,
                                                std::vector<double> weights, double scaling_factor,
                                                const std::vector<bool> &bin_mask1,
                                                const std::vector<bool> &bin_mask2,
                                                std::uint64_t min_delta_, std::uint64_t max_delta_)
    : _chrom1(std::move(chrom1)),
      _chrom2(std::move(chrom2)),
      _bins(std::move(bins)),
      _weights(std::move(weights)),
      _min_delta(min_delta_),
      _max_delta(max_delta_) {
  if (scaling_factor != 1) {
    std::ranges::transform(_weights, _weights.begin(), [&](const auto n) {
      const auto scaled_n = n / scaling_factor;
      return std::isfinite(scaled_n) ? scaled_n : 0.0;
    });
  }

  MatrixStats<double> stats(_chrom1, _chrom2, bin_mask1, bin_mask2, _bins.resolution(), min_delta_,
                            max_delta_, _weights);
  stats.add(pixels);

  _marginals1 = std::move(stats.marginals1);
  _marginals2 = std::move(stats.marginals2);
  _nnz = stats.nnz;
  _sum = stats.sum;

  assert(internal::is_close(_sum, std::accumulate(_marginals1->begin(), _marginals1->end(), 0.0)));
  assert(internal::is_close(_sum, std::accumulate(_marginals2->begin(), _marginals2->end(), 0.0)));
}

template <typename Pixels, typename PixelsGW>
  requires PixelRange<Pixels> && PixelRange<PixelsGW>
inline ExpectedMatrixStats::ExpectedMatrixStats(const Pixels &pixels, const PixelsGW &pixels_gw,
                                                const hictk::Chromosome &chrom1,
                                                const hictk::Chromosome &chrom2,
                                                const hictk::BinTable &bins,
                                                const std::vector<bool> &bin_mask1,
                                                const std::vector<bool> &bin_mask2,
                                                std::uint64_t min_delta_, std::uint64_t max_delta_)
    : ExpectedMatrixStats(pixels, chrom1, chrom2, bins,
                          compute_weights(pixels_gw, chrom1, chrom2, bins, min_delta_, max_delta_),
                          1.0, bin_mask1, bin_mask2, min_delta_, max_delta_) {}

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
ExpectedMatrixStats::build_expected_vector(const Pixels &pixels, const hictk::BinTable &bins,
                                           std::uint64_t min_delta_, std::uint64_t max_delta_) {
  if (std::ranges::empty(pixels)) {
    phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
    std::ranges::transform(
        bins.chromosomes(), std::inserter(scaling_factors, scaling_factors.begin()),
        [](const hictk::Chromosome &chrom) { return std::make_pair(chrom, 0.0); });
    return {std::vector<double>(bins.size(), 0), scaling_factors};
  }

  const auto bins_ = std::make_shared<const hictk::BinTable>(bins);

  hictk::ExpectedValuesAggregator aggr(bins_);

  for (const auto &p : pixels) {
    const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
    if (delta >= min_delta_ && delta < max_delta_) [[likely]] {
      aggr.add(p);
    }
  }

  aggr.compute_density();

  auto weights = aggr.weights();
  const auto &chrom = bins.chromosomes().longest_chromosome();
  weights.resize((chrom.size() + bins.resolution() - 1) / bins.resolution(),
                 std::numeric_limits<double>::quiet_NaN());

  return {weights, aggr.scaling_factors()};
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::vector<double> ExpectedMatrixStats::compute_weights(
    const Pixels &pixels, const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
    const hictk::BinTable &bins, std::uint64_t min_delta_, std::uint64_t max_delta_) {
  if (chrom1 != chrom2) {
    return {};
  }

  auto [weights, scaling_factors] = build_expected_vector(pixels, bins, min_delta_, max_delta_);

  weights.resize((chrom1.size() + bins.resolution() - 1) / bins.resolution(),
                 std::numeric_limits<double>::quiet_NaN());

  const auto sf = scaling_factors.at(chrom1);
  std::ranges::transform(weights, weights.begin(), [&](const auto n) { return n / sf; });

  return weights;
}

}  // namespace nchg
