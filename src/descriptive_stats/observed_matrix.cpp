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

#include "nchg/observed_matrix.hpp"

#include <cstddef>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <vector>

namespace nchg {

std::uint32_t ObservedMatrix::resolution() const noexcept { return _bins.resolution(); }

std::size_t ObservedMatrix::num_rows() const noexcept { return _bins.subset(chrom1()).size(); }

std::size_t ObservedMatrix::num_cols() const noexcept { return _bins.subset(chrom2()).size(); }

const hictk::Chromosome &ObservedMatrix::chrom1() const noexcept { return _chrom1; }

const hictk::Chromosome &ObservedMatrix::chrom2() const noexcept { return _chrom2; }

std::uint64_t ObservedMatrix::nnz() const noexcept { return _nnz; }

std::uint64_t ObservedMatrix::sum() const noexcept { return _sum; }

double ObservedMatrix::nnz_avg() const noexcept {
  if (nnz() == 0) [[unlikely]] {
    return 0.0;
  }
  return static_cast<double>(sum()) / static_cast<double>(nnz());
}

double ObservedMatrix::mad_max() const noexcept { return _mad_max; }

std::uint64_t ObservedMatrix::min_delta() const noexcept { return _min_delta; }

std::uint64_t ObservedMatrix::max_delta() const noexcept { return _max_delta; }

const std::vector<std::uint64_t> &ObservedMatrix::marginals1() const noexcept {
  return *_marginals1;
}

const std::vector<std::uint64_t> &ObservedMatrix::marginals2() const noexcept {
  return *_marginals2;
}

}  // namespace nchg
