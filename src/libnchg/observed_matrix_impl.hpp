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

#include <cassert>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <memory>
#include <utility>
#include <vector>

#include "nchg/concepts.hpp"

namespace nchg {

template <typename Pixels>
  requires PixelRange<Pixels>
inline auto ObservedMatrix::compute_stats(const Pixels &pixels, const hictk::Chromosome &chrom1,
                                          const hictk::Chromosome &chrom2,
                                          const hictk::BinTable &bins,
                                          const std::vector<bool> &bin_mask1,
                                          const std::vector<bool> &bin_mask2,
                                          std::uint64_t min_delta_, std::uint64_t max_delta_) {
  assert(min_delta_ <= max_delta_);
  struct Result {
    std::shared_ptr<std::vector<std::uint64_t>> marginals1{
        std::make_shared<std::vector<std::uint64_t>>()};
    std::shared_ptr<std::vector<std::uint64_t>> marginals2{
        std::make_shared<std::vector<std::uint64_t>>()};
    std::uint64_t sum{};
    std::uint64_t nnz{};
  };

  Result res{};
  const auto intra_matrix = chrom1 == chrom2;
  if (intra_matrix) {
    res.marginals2 = res.marginals1;
  }

  const auto bins1 = bins.subset(chrom1);
  const auto bins2 = bins.subset(chrom2);

  if (bins1.empty() || bins2.empty()) {
    return res;
  }

  const auto num_bins1 = bins1.size();
  const auto num_bins2 = bins2.size();

  auto &marginals1_ = *res.marginals1;
  auto &marginals2_ = *res.marginals2;
  auto &sum_ = res.sum;
  auto &nnz_ = res.nnz;
  marginals1_.resize(num_bins1, 0);
  marginals2_.resize(num_bins2, 0);

  for (const hictk::Pixel<N> &p : pixels) {
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

    if (intra_matrix) {
      sum_ += bin1 == bin2 ? p.count : 2 * p.count;
      nnz_ += bin1 == bin2 ? 1ULL : 2ULL;
    } else {
      sum_ += p.count;
      ++nnz_;
    }

    const auto i1 = bin1.rel_id();
    marginals1_[i1] += p.count;

    if (!intra_matrix || bin1 != bin2) {
      const auto i2 = bin2.rel_id();
      marginals2_[i2] += p.count;
    }
  };

  return res;
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline ObservedMatrix::ObservedMatrix(const Pixels &pixels, hictk::Chromosome chrom1,
                                      hictk::Chromosome chrom2, hictk::BinTable bins,
                                      double mad_max_, const std::vector<bool> &bin_mask1,
                                      const std::vector<bool> &bin_mask2, std::uint64_t min_delta_,
                                      std::uint64_t max_delta_)
    : _chrom1(std::move(chrom1)),
      _chrom2(std::move(chrom2)),
      _bins(std::move(bins)),
      _mad_max(mad_max_),
      _min_delta(min_delta_),
      _max_delta(max_delta_) {
  auto stats =
      compute_stats(pixels, _chrom1, _chrom2, _bins, bin_mask1, bin_mask2, min_delta_, max_delta_);
  _marginals1 = std::move(stats.marginals1);
  _marginals2 = std::move(stats.marginals2);
  _nnz = stats.nnz;
  _sum = stats.sum;
}

}  // namespace nchg
