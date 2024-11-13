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

#include <cstddef>
#include <cstdint>
#include <functional>
#include <hictk/bin.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>

#include "nchg/concepts.hpp"

namespace nchg {

template <typename N>
  requires arithmetic<N>
struct MatrixStats {
  using T = std::conditional_t<std::is_floating_point_v<N>, N, std::uint64_t>;
  using MarginalBuff = std::vector<T>;

  std::shared_ptr<MarginalBuff> marginals1{};
  std::shared_ptr<MarginalBuff> marginals2{};
  T sum{};
  std::uint64_t nnz{};

 private:
  bool _intra_matrix{};
  std::span<const double> _weights{};
  const std::vector<bool> *_bin_mask1{};
  const std::vector<bool> *_bin_mask2{};
  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{std::numeric_limits<std::uint64_t>::max()};

 public:
  MatrixStats() = default;
  MatrixStats(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
              const std::vector<bool> &bin_mask1, const std::vector<bool> &bin_mask2,
              std::uint32_t bin_size, std::uint64_t min_delta, std::uint64_t max_delta,
              std::span<const double> weights = {});

  void add(const hictk::Pixel<N> &p);

  template <typename Pixels>
    requires PixelRange<Pixels>
  void add(const Pixels &pixels);

 private:
  [[nodiscard]] static constexpr std::size_t compute_num_bins(const hictk::Chromosome &chrom,
                                                              std::uint32_t bin_size) noexcept;
  [[nodiscard]] bool is_masked(const hictk::Bin &bin1, const hictk::Bin &bin2) const noexcept;
  constexpr auto get_count(const hictk::Pixel<N> &p) const noexcept -> T;
};

}  // namespace nchg

#include "../../matrix_stats_impl.hpp"
