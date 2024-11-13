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
#include "nchg/matrix_stats.hpp"

namespace nchg {

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
  MatrixStats<N> stats(_chrom1, _chrom2, bin_mask1, bin_mask2, bins.resolution(), _min_delta,
                       _max_delta);
  stats.add(pixels);

  _marginals1 = std::move(stats.marginals1);
  _marginals2 = std::move(stats.marginals2);
  _nnz = stats.nnz;
  _sum = stats.sum;

  assert(_nnz <= _sum);
  assert(_sum == std::accumulate(_marginals1->begin(), _marginals1->end(), std::uint64_t{}));
  assert(_sum == std::accumulate(_marginals2->begin(), _marginals2->end(), std::uint64_t{}));
}

}  // namespace nchg
