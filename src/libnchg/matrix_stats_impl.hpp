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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <hictk/bin.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/concepts.hpp"

namespace nchg {

template <typename N>
  requires arithmetic<N>
MatrixStats<N>::MatrixStats(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
                            const std::vector<bool> &bin_mask1, const std::vector<bool> &bin_mask2,
                            std::uint32_t bin_size, std::uint64_t min_delta,
                            std::uint64_t max_delta, std::span<const double> weights)
    : marginals1(std::make_shared<MarginalBuff>(compute_num_bins(chrom1, bin_size), 0)),
      marginals2(chrom1 == chrom2
                     ? marginals1
                     : std::make_shared<MarginalBuff>(compute_num_bins(chrom2, bin_size), 0)),
      _intra_matrix(chrom1 == chrom2),
      _weights(weights),
      _bin_mask1(&bin_mask1),
      _bin_mask2(&bin_mask2),
      _min_delta(min_delta),
      _max_delta(max_delta) {
  if (_min_delta > _max_delta) {
    throw std::logic_error("min_delta must be smaller or equal to max delta");
  }
}

template <typename N>
  requires arithmetic<N>
[[nodiscard]] bool MatrixStats<N>::is_masked(const hictk::Bin &bin1,
                                             const hictk::Bin &bin2) const noexcept {
  assert(bin1 <= bin2);
  assert(_bin_mask1);
  assert(_bin_mask2);

  const auto delta = bin2.start() - bin1.start();
  if (_intra_matrix && (delta < _min_delta || delta >= _max_delta)) [[unlikely]] {
    return true;
  }

  const auto bin1_id = bin1.rel_id();
  const auto bin2_id = bin2.rel_id();

  assert(bin1_id < _bin_mask1->size());
  assert(bin2_id < _bin_mask2->size());

  const auto bin1_masked = !_bin_mask1->empty() && (*_bin_mask1)[bin1_id];
  const auto bin2_masked = !_bin_mask2->empty() && (*_bin_mask2)[bin2_id];

  if (bin1_masked || bin2_masked) [[unlikely]] {
    return true;
  }

  return false;
}

template <typename N>
  requires arithmetic<N>
constexpr std::size_t MatrixStats<N>::compute_num_bins(const hictk::Chromosome &chrom,
                                                       std::uint32_t bin_size) noexcept {
  assert(bin_size != 0);
  return (chrom.size() + bin_size - 1) / bin_size;
}

template <typename N>
  requires arithmetic<N>
constexpr auto MatrixStats<N>::get_count(const hictk::Pixel<N> &p) const noexcept -> T {
  const auto &[bin1, bin2] = p.coords;

  const auto count = [&] {
    if constexpr (std::is_floating_point_v<T>) {
      if (_intra_matrix && !_weights.empty()) {
        assert(bin1 <= bin2);
        const auto diag = bin2.id() - bin1.id();
        assert(diag < _weights.size());
        return _weights[diag];
      }
    }

    return conditional_static_cast<T>(p.count);
  }();

  if constexpr (std::is_floating_point_v<T>) {
    if (std::isnan(count)) [[unlikely]] {
      return T{0};
    }
  }

  return count;
}

template <typename N>
  requires arithmetic<N>
inline void MatrixStats<N>::add(const hictk::Pixel<N> &p) {
  const auto &[bin1, bin2] = p.coords;

  if (is_masked(bin1, bin2)) [[unlikely]] {
    return;
  }

  const auto count = get_count(p);

  if (_intra_matrix) {
    sum += bin1 == bin2 ? count : 2 * count;
    nnz += bin1 == bin2 ? 1ULL : 2ULL;
  } else {
    sum += count;
    ++nnz;
  }

  const auto i1 = bin1.rel_id();
  (*marginals1)[i1] += count;

  if (!_intra_matrix || bin1 != bin2) {
    const auto i2 = bin2.rel_id();
    (*marginals2)[i2] += count;
  }
}

template <typename N>
  requires arithmetic<N>
template <typename Pixels>
  requires PixelRange<Pixels>
inline void MatrixStats<N>::add(const Pixels &pixels) {
  for (const auto &p : pixels) {
    add(p);
  }
}

}  // namespace nchg
