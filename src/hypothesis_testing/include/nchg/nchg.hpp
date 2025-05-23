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
#include <parallel_hashmap/phmap.h>

#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/cooler/pixel_selector.hpp>
#include <hictk/file.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/hic/pixel_selector.hpp>
#include <hictk/pixel.hpp>
#include <memory>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/genomic_domains.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {

struct NCHGResult {
  hictk::Pixel<std::uint64_t> pixel;
  double expected{};
  double pval{};
  double log_ratio{};
  double odds_ratio{};
  double omega{};

  bool operator<(const NCHGResult& other) const noexcept;
  bool operator>(const NCHGResult& other) const noexcept;
  bool operator==(const NCHGResult& other) const noexcept;
  bool operator!=(const NCHGResult& other) const noexcept;
};

class NCHG {
 public:
  using N = std::uint32_t;
  template <typename PixelIt>
  class iterator;
  using Stats = NCHGResult;

 private:
  std::shared_ptr<const hictk::File> _fp;

  hictk::Chromosome _chrom1;
  hictk::Chromosome _chrom2;

  std::shared_ptr<const ExpectedMatrixStats> _exp_matrix;
  std::shared_ptr<const ObservedMatrix> _obs_matrix;
  ExpectedValues _expected_values;

  mutable std::vector<double> _nchg_pval_buffer;

 public:
  using IteratorVariant =
      std::variant<iterator<hictk::cooler::PixelSelector>, iterator<hictk::hic::PixelSelector>,
                   iterator<hictk::hic::PixelSelectorAll>>;

  using Params = ExpectedValues::Params;
  static constexpr auto DefaultParams = ExpectedValues::DefaultParams;
  NCHG(const std::shared_ptr<const hictk::File>& f, const hictk::Chromosome& chrom1,
       const hictk::Chromosome& chrom2, const Params& params);
  NCHG(std::shared_ptr<const hictk::File> f, hictk::Chromosome chrom1, hictk::Chromosome chrom2,
       ExpectedValues expected_values);

  [[nodiscard]] auto params() const noexcept -> Params;

  [[nodiscard]] const ObservedMatrix& observed_matrix() const noexcept;
  [[nodiscard]] const ExpectedMatrixStats& expected_matrix() const noexcept;

  [[nodiscard]] auto compute(const BEDPE& domain, std::uint64_t obs, double exp,
                             double bad_bin_fraction) const -> Stats;

  [[nodiscard]] auto cbegin(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> IteratorVariant;
  [[nodiscard]] auto cend(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> IteratorVariant;

  [[nodiscard]] auto begin(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> IteratorVariant;
  [[nodiscard]] auto end(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> IteratorVariant;

 private:
  [[nodiscard]] static auto init_matrices(const hictk::Chromosome& chrom1,
                                          const hictk::Chromosome& chrom2, const hictk::File& f,
                                          const ExpectedValues& expected_values,
                                          const std::vector<bool>& bin1_mask,
                                          const std::vector<bool>& bin2_mask, double mad_max_,
                                          std::uint64_t min_delta_, std::uint64_t max_delta_)
      -> std::pair<std::shared_ptr<const ObservedMatrix>,
                   std::shared_ptr<const ExpectedMatrixStats>>;

  [[nodiscard]] static double compute_cumulative_nchg(std::vector<double>& buffer,
                                                      std::uint64_t obs, std::uint64_t N1,
                                                      std::uint64_t N2, std::uint64_t N,
                                                      double odds, double precision,
                                                      bool lower_tail = false);

  [[nodiscard]] static double compute_pvalue_nchg(std::vector<double>& buffer, std::uint64_t obs,
                                                  std::uint64_t N1, std::uint64_t N2,
                                                  std::uint64_t N, double odds,
                                                  double precision = 1.0e-20,
                                                  double min_omega = 0.1);

  [[nodiscard]] static double compute_pvalue_nchg(std::uint64_t obs, std::uint64_t N1,
                                                  std::uint64_t N2, std::uint64_t N, double odds,
                                                  double precision = 1.0e-20,
                                                  double min_omega = 0.1);

  [[nodiscard]] std::uint64_t compute_N1(const hictk::GenomicInterval& range1,
                                         const hictk::GenomicInterval& range2,
                                         double max_bad_bin_threshold) const noexcept;
  [[nodiscard]] std::uint64_t compute_N2(const hictk::GenomicInterval& range1,
                                         const hictk::GenomicInterval& range2,
                                         double max_bad_bin_threshold) const noexcept;
  [[nodiscard]] double compute_L1(const hictk::GenomicInterval& range1,
                                  const hictk::GenomicInterval& range2,
                                  double max_bad_bin_threshold) const noexcept;
  [[nodiscard]] double compute_L2(const hictk::GenomicInterval& range1,
                                  const hictk::GenomicInterval& range2,
                                  double max_bad_bin_threshold) const noexcept;
  [[nodiscard]] auto aggregate_pixels(const hictk::GenomicInterval& range1,
                                      const hictk::GenomicInterval& range2) const;

  [[nodiscard]] static NCHGResult compute_stats(hictk::Pixel<std::uint64_t> pixel, double exp,
                                                std::uint64_t obs_sum, std::uint64_t N1,
                                                std::uint64_t N2, double exp_sum, double L1,
                                                double L2, std::vector<double>& buffer);

 public:
  template <typename PixelSelector>
  class iterator {
    using PixelIt = decltype(std::declval<PixelSelector>().template begin<N>());
    std::shared_ptr<const PixelSelector> _sel;
    PixelIt _pixel_it{};
    PixelIt _sentinel_it{};

    std::shared_ptr<const ObservedMatrix> _obs;
    std::shared_ptr<const ExpectedMatrixStats> _exp;

    std::shared_ptr<const std::vector<bool>> _bin_mask1;
    std::shared_ptr<const std::vector<bool>> _bin_mask2;

    std::shared_ptr<std::vector<double>> _buffer{std::make_shared<std::vector<double>>()};

    std::uint64_t _min_delta{};
    std::uint64_t _max_delta{};

    mutable Stats _value{};
    mutable bool _read_value{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Stats;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelSelector sel, std::shared_ptr<const ObservedMatrix> obs,
             std::shared_ptr<const ExpectedMatrixStats> exp,
             std::shared_ptr<const std::vector<bool>> bin_mask1,
             std::shared_ptr<const std::vector<bool>> bin_mask2, std::uint64_t min_delta,
             std::uint64_t max_delta);

    [[nodiscard]] static auto at_end(PixelSelector sel, std::shared_ptr<const ObservedMatrix> obs,
                                     std::shared_ptr<const ExpectedMatrixStats> exp) -> iterator;

    iterator(const iterator& other);
    iterator(iterator&& other) noexcept = default;

    ~iterator() noexcept = default;

    auto operator=(const iterator& other) -> iterator&;
    auto operator=(iterator&& other) noexcept -> iterator& = default;

    [[nodiscard]] bool operator==(const iterator& other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator& other) const noexcept;

    [[nodiscard]] bool operator<(const iterator& other) const noexcept;
    [[nodiscard]] bool operator<=(const iterator& other) const noexcept;

    [[nodiscard]] bool operator>(const iterator& other) const noexcept;
    [[nodiscard]] bool operator>=(const iterator& other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator&;
    auto operator++(int) -> iterator;

   private:
    void jump_to_next_valid_pixel();
    [[nodiscard]] constexpr std::uint64_t compute_N1(
        const hictk::Pixel<std::uint64_t>& pixel) const noexcept;
    [[nodiscard]] constexpr std::uint64_t compute_N2(
        const hictk::Pixel<std::uint64_t>& pixel) const noexcept;
    [[nodiscard]] constexpr double compute_L1(
        const hictk::Pixel<std::uint64_t>& pixel) const noexcept;
    [[nodiscard]] constexpr double compute_L2(
        const hictk::Pixel<std::uint64_t>& pixel) const noexcept;
  };
};

}  // namespace nchg

#include "../../nchg_impl.hpp"
