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

#include "nchg/expected_matrix.hpp"

#include <parallel_hashmap/btree.h>

#include <cstddef>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <vector>

namespace nchg {

std::uint32_t ExpectedMatrixStats::resolution() const noexcept { return _bins.resolution(); }

std::size_t ExpectedMatrixStats::num_rows() const noexcept { return _bins.subset(chrom1()).size(); }

std::size_t ExpectedMatrixStats::num_cols() const noexcept { return _bins.subset(chrom2()).size(); }

const hictk::Chromosome &ExpectedMatrixStats::chrom1() const noexcept { return _chrom1; }

const hictk::Chromosome &ExpectedMatrixStats::chrom2() const noexcept { return _chrom2; }

std::uint64_t ExpectedMatrixStats::nnz() const noexcept { return _nnz; }

double ExpectedMatrixStats::sum() const noexcept { return _sum; }

double ExpectedMatrixStats::nnz_avg() const noexcept {
  if (nnz() == 0) [[unlikely]] {
    return 0.0;
  }
  return sum() / static_cast<double>(nnz());
}

const std::vector<double> &ExpectedMatrixStats::weights() const noexcept { return _weights; }

const phmap::btree_map<hictk::Chromosome, double> &ExpectedMatrixStats::scaling_factors()
    const noexcept {
  return _scaling_factors;
}

std::uint64_t ExpectedMatrixStats::min_delta() const noexcept { return _min_delta; }

std::uint64_t ExpectedMatrixStats::max_delta() const noexcept { return _max_delta; }

double ExpectedMatrixStats::at(std::uint64_t i, std::uint64_t j) const {
  if (chrom1() == chrom2()) {
    return _weights.at(j - i);
  }
  return nnz_avg();
}

const std::vector<double> &ExpectedMatrixStats::marginals1() const noexcept { return *_marginals1; }
const std::vector<double> &ExpectedMatrixStats::marginals2() const noexcept { return *_marginals2; }

}  // namespace nchg
