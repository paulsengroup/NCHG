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

#include <cstddef>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>
#include <limits>
#include <memory>
#include <vector>

#include "nchg/concepts.hpp"

namespace nchg {

class ExpectedMatrixStats {
  hictk::Chromosome _chrom1{};
  hictk::Chromosome _chrom2{};

  hictk::BinTable _bins{};
  std::vector<double> _weights{};
  phmap::btree_map<hictk::Chromosome, double> _scaling_factors{};

  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{};

  std::shared_ptr<std::vector<double>> _marginals1{};
  std::shared_ptr<std::vector<double>> _marginals2{};

  double _sum{};
  std::uint64_t _nnz{};

 public:
  template <typename Pixels>
    requires PixelRange<Pixels>
  ExpectedMatrixStats(const Pixels &pixels, hictk::Chromosome chrom1, hictk::Chromosome chrom2,
                      hictk::BinTable bins, std::vector<double> weights, double scaling_factor,
                      const std::vector<bool> &bin_mask1 = {},
                      const std::vector<bool> &bin_mask2 = {}, std::uint64_t min_delta_ = 0,
                      std::uint64_t max_delta_ = std::numeric_limits<std::uint64_t>::max());

  template <typename Pixels, typename PixelsGW>
    requires PixelRange<Pixels> && PixelRange<PixelsGW>
  ExpectedMatrixStats(const Pixels &pixels, const PixelsGW &pixels_gw,
                      const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
                      const hictk::BinTable &bins, const std::vector<bool> &bin_mask1 = {},
                      const std::vector<bool> &bin_mask2 = {}, std::uint64_t min_delta_ = 0,
                      std::uint64_t max_delta_ = std::numeric_limits<std::uint64_t>::max());

  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] std::size_t num_rows() const noexcept;
  [[nodiscard]] std::size_t num_cols() const noexcept;

  [[nodiscard]] const hictk::Chromosome &chrom1() const noexcept;
  [[nodiscard]] const hictk::Chromosome &chrom2() const noexcept;

  [[nodiscard]] std::uint64_t nnz() const noexcept;
  [[nodiscard]] double sum() const noexcept;
  [[nodiscard]] double nnz_avg() const noexcept;

  [[nodiscard]] const std::vector<double> &weights() const noexcept;
  [[nodiscard]] const phmap::btree_map<hictk::Chromosome, double> &scaling_factors() const noexcept;

  [[nodiscard]] std::uint64_t min_delta() const noexcept;
  [[nodiscard]] std::uint64_t max_delta() const noexcept;

  [[nodiscard]] double at(std::uint64_t i, std::uint64_t j) const;

  [[nodiscard]] const std::vector<double> &marginals1() const noexcept;
  [[nodiscard]] const std::vector<double> &marginals2() const noexcept;

  template <typename Pixels>
    requires PixelRange<Pixels>
  [[nodiscard]] static std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
  build_expected_vector(const Pixels &pixels, const hictk::BinTable &bins, std::uint64_t min_delta_,
                        std::uint64_t max_delta_);

  template <typename Pixels>
    requires PixelRange<Pixels>
  static std::vector<double> compute_weights(const Pixels &pixels, const hictk::Chromosome &chrom1,
                                             const hictk::Chromosome &chrom2,
                                             const hictk::BinTable &bins, std::uint64_t min_delta_,
                                             std::uint64_t max_delta_);
};

}  // namespace nchg

#include "../../expected_matrix_impl.hpp"
