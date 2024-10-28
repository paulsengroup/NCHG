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

#include <algorithm>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/expected_values_aggregator.hpp>
#include <hictk/file.hpp>
#include <hictk/pixel.hpp>
#include <utility>
#include <variant>
#include <vector>

#include "nchg/concepts.hpp"

namespace nchg {

template <typename Pixels>
  requires PixelRange<Pixels>
inline auto ExpectedMatrix::compute_stats(const Pixels &pixels, const hictk::Chromosome &chrom1,
                                          const hictk::Chromosome &chrom2,
                                          const hictk::BinTable &bins,
                                          const std::vector<double> &weights,
                                          const std::vector<bool> &bin_mask1,
                                          const std::vector<bool> &bin_mask2,
                                          std::uint64_t min_delta_, std::uint64_t max_delta_) {
  struct Result {
    std::shared_ptr<std::vector<double>> marginals1{std::make_shared<std::vector<double>>()};
    std::shared_ptr<std::vector<double>> marginals2{std::make_shared<std::vector<double>>()};
    double sum{};
    std::uint64_t nnz{};
  };

  Result res{};
  const auto intra_matrix = chrom1 == chrom2;
  if (intra_matrix) {
    res.marginals2 = res.marginals1;
  }

  const auto bin_size = bins.resolution();
  const auto num_bins1 = (chrom1.size() + bin_size - 1) / bin_size;
  const auto num_bins2 = (chrom2.size() + bin_size - 1) / bin_size;

  if (num_bins1 == 0 || num_bins2 == 0) {
    return res;
  }

  auto &marginals1_ = *res.marginals1;
  auto &marginals2_ = *res.marginals2;
  auto &sum_ = res.sum;
  auto &nnz_ = res.nnz;
  marginals1_.resize(num_bins1, 0);
  marginals2_.resize(num_bins2, 0);

  for (const auto &p : pixels) {
    const auto &bin1 = p.coords.bin1;
    const auto &bin2 = p.coords.bin2;
    const auto delta = bin2.start() - bin1.start();
    if (intra_matrix && (delta < min_delta_ || delta >= max_delta_)) {
      continue;
    }

    const auto bin1_id = p.coords.bin1.rel_id();
    const auto bin2_id = p.coords.bin2.rel_id();

    const auto bin1_masked = !bin_mask1.empty() && bin_mask1[bin1_id];
    const auto bin2_masked = !bin_mask2.empty() && bin_mask2[bin2_id];

    if (bin1_masked || bin2_masked) {
      continue;
    }

    auto count =
        intra_matrix ? weights.at(bin2.id() - bin1.id()) : conditional_static_cast<double>(p.count);
    if (std::isnan(count)) {
      count = 0.0;
    }

    if (intra_matrix) {
      sum_ += bin1 == bin2 ? count : 2 * count;
      nnz_ += bin1 == bin2 ? 1ULL : 2ULL;
    } else {
      sum_ += count;
      nnz_++;
    }

    const auto i1 = bin1.rel_id();
    marginals1_[i1] += count;

    if (!intra_matrix || bin1 != bin2) {
      const auto i2 = bin2.rel_id();
      marginals2_[i2] += count;
    }
  }

  return res;
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline ExpectedMatrix::ExpectedMatrix(const Pixels &pixels, hictk::Chromosome chrom1,
                                      hictk::Chromosome chrom2, hictk::BinTable bins,
                                      std::vector<double> weights, double scaling_factor,
                                      const std::vector<bool> &bin_mask1,
                                      const std::vector<bool> &bin_mask2, std::uint64_t min_delta_,
                                      std::uint64_t max_delta_)
    : _chrom1(std::move(chrom1)),
      _chrom2(std::move(chrom2)),
      _bins(std::move(bins)),
      _weights(std::move(weights)),
      _min_delta(min_delta_),
      _max_delta(max_delta_) {
  if (scaling_factor != 1) {
    std::transform(_weights.begin(), _weights.end(), _weights.begin(),
                   [&](const auto n) { return n / scaling_factor; });
  }

  auto stats = compute_stats(pixels, _chrom1, _chrom2, _bins, _weights, bin_mask1, bin_mask2,
                             min_delta_, max_delta_);
  _marginals1 = std::move(stats.marginals1);
  _marginals2 = std::move(stats.marginals2);
  _nnz = stats.nnz;
  _sum = stats.sum;
}

inline std::uint32_t ExpectedMatrix::resolution() const noexcept { return _bins.resolution(); }

inline std::size_t ExpectedMatrix::num_rows() const noexcept {
  return _bins.subset(chrom1()).size();
}

inline std::size_t ExpectedMatrix::num_cols() const noexcept {
  return _bins.subset(chrom2()).size();
}

inline const hictk::Chromosome &ExpectedMatrix::chrom1() const noexcept { return _chrom1; }

inline const hictk::Chromosome &ExpectedMatrix::chrom2() const noexcept { return _chrom2; }

inline std::uint64_t ExpectedMatrix::nnz() const noexcept { return _nnz; }

inline double ExpectedMatrix::sum() const noexcept { return _sum; }

inline double ExpectedMatrix::nnz_avg() const noexcept {
  return sum() / static_cast<double>(nnz());
}

inline const std::vector<double> &ExpectedMatrix::weights() const noexcept { return _weights; }

inline const phmap::btree_map<hictk::Chromosome, double> &ExpectedMatrix::scaling_factors()
    const noexcept {
  return _scaling_factors;
}

inline std::uint64_t ExpectedMatrix::min_delta() const noexcept { return _min_delta; }

inline std::uint64_t ExpectedMatrix::max_delta() const noexcept { return _max_delta; }

inline double ExpectedMatrix::at(std::uint64_t i, std::uint64_t j) const {
  if (chrom1() == chrom2()) {
    return _weights.at(j - i);
  }
  return nnz_avg();
}

inline const std::vector<double> &ExpectedMatrix::marginals1() const noexcept {
  return *_marginals1;
}
inline const std::vector<double> &ExpectedMatrix::marginals2() const noexcept {
  return *_marginals2;
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
ExpectedMatrix::build_expected_vector(const Pixels &pixels, const hictk::BinTable &bins,
                                      std::uint64_t min_delta_, std::uint64_t max_delta_) {
  if (std::ranges::empty(pixels)) {
    phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
    std::transform(bins.chromosomes().begin(), bins.chromosomes().end(),
                   std::inserter(scaling_factors, scaling_factors.begin()),
                   [](const hictk::Chromosome &chrom) { return std::make_pair(chrom, 0.0); });
    return std::make_pair(std::vector<double>(bins.size(), 0), scaling_factors);
  }

  const auto bins_ = std::make_shared<const hictk::BinTable>(bins);

  hictk::ExpectedValuesAggregator aggr(bins_);

  for (const auto &p : pixels) {
    const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
    if (delta >= min_delta_ && delta < max_delta_) {
      aggr.add(p);
    }
  }

  aggr.compute_density();

  auto weights = aggr.weights();
  const auto &chrom = bins.chromosomes().longest_chromosome();
  weights.resize((chrom.size() + bins.resolution() - 1) / bins.resolution(),
                 std::numeric_limits<double>::quiet_NaN());

  return std::make_pair(std::move(weights), aggr.scaling_factors());
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::vector<double> ExpectedMatrix::compute_weights(
    const Pixels &pixels, const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
    const hictk::BinTable &bins, std::uint64_t min_delta_, std::uint64_t max_delta_) {
  if (chrom1 != chrom2) {
    return {};
  }

  auto [weights, scaling_factors] = build_expected_vector(pixels, bins, min_delta_, max_delta_);

  weights.resize((chrom1.size() + bins.resolution() - 1) / bins.resolution(),
                 std::numeric_limits<double>::quiet_NaN());

  const auto sf = scaling_factors.at(chrom1);
  std::ranges::transform(weights, weights, [&](const auto n) { return n / sf; });

  return weights;
}

}  // namespace nchg
