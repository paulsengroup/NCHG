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

#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>
#include <limits>
#include <memory>
#include <vector>

#include "nchg/concepts.hpp"

namespace nchg {

class ObservedMatrix {
  using N = std::uint32_t;
  using MarginalBuff = std::vector<std::uint64_t>;

  hictk::Chromosome _chrom1{};
  hictk::Chromosome _chrom2{};

  hictk::BinTable _bins{};

  std::shared_ptr<const MarginalBuff> _marginals1{};
  std::shared_ptr<const MarginalBuff> _marginals2{};

  double _mad_max{};
  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{};

  std::uint64_t _nnz{};
  std::uint64_t _sum{};

 public:
  ObservedMatrix() = delete;

  template <typename Pixels>
    requires PixelRange<Pixels>
  ObservedMatrix(const Pixels& pixels, hictk::Chromosome chrom1, hictk::Chromosome chrom2,
                 hictk::BinTable bins, double mad_max_ = 0.0,
                 const std::vector<bool>& bin_mask1 = {}, const std::vector<bool>& bin_mask2 = {},
                 std::uint64_t min_delta_ = 0,
                 std::uint64_t max_delta_ = std::numeric_limits<std::uint64_t>::max());

  ObservedMatrix(hictk::Chromosome chrom1, hictk::Chromosome chrom2, hictk::BinTable bins,
                 std::shared_ptr<const MarginalBuff> marginals1_,
                 std::shared_ptr<const MarginalBuff> marginals2_, std::uint64_t nnz_,
                 std::uint64_t sum_, double mad_max_ = 0.0, std::uint64_t min_delta_ = 0,
                 std::uint64_t max_delta_ = std::numeric_limits<std::uint64_t>::max()) noexcept;

  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] std::size_t num_rows() const noexcept;
  [[nodiscard]] std::size_t num_cols() const noexcept;

  [[nodiscard]] const hictk::Chromosome& chrom1() const noexcept;
  [[nodiscard]] const hictk::Chromosome& chrom2() const noexcept;

  [[nodiscard]] std::uint64_t nnz() const noexcept;
  [[nodiscard]] std::uint64_t sum() const noexcept;
  [[nodiscard]] double nnz_avg() const noexcept;

  [[nodiscard]] double mad_max() const noexcept;
  [[nodiscard]] std::uint64_t min_delta() const noexcept;
  [[nodiscard]] std::uint64_t max_delta() const noexcept;

  [[nodiscard]] const std::vector<std::uint64_t>& marginals1() const noexcept;
  [[nodiscard]] const std::vector<std::uint64_t>& marginals2() const noexcept;

 private:
  template <typename Pixels>
    requires PixelRange<Pixels>
  static auto compute_stats(const Pixels& pixels, const hictk::Chromosome& chrom1,
                            const hictk::Chromosome& chrom2, const hictk::BinTable& bins,
                            const std::vector<bool>& bin_mask1, const std::vector<bool>& bin_mask2,
                            std::uint64_t min_delta_, std::uint64_t max_delta_);
};

}  // namespace nchg

#include "../../observed_matrix_impl.hpp"
