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
#include <hictk/cooler/cooler.hpp>
#include <hictk/file.hpp>
#include <hictk/hic.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <memory>

#include "nchg/expected_matrix.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {

template <typename File>
class NCHG {
 public:
  class iterator;

  struct Stats {
    hictk::Pixel<std::uint32_t> pixel{};
    double expected{};
    double pval{};
    double odds_ratio{};
    double omega{};
  };

 private:
  using N = std::uint32_t;
  using ThinPixelIt = decltype(std::declval<File>().fetch("chr1", "chr2").template begin<N>());
  using PixelIt =
      decltype(std::declval<hictk::transformers::JoinGenomicCoords<ThinPixelIt>>().begin());

  std::shared_ptr<const File> _fp;

  hictk::Chromosome _chrom1{};
  hictk::Chromosome _chrom2{};

  std::shared_ptr<const ExpectedMatrix<PixelIt>> _exp_matrix{};
  std::shared_ptr<const ObservedMatrix<PixelIt>> _obs_matrix{};
  ExpectedValues<File> _expected_values{};

  mutable std::vector<double> _nchg_pval_buffer{};

 public:
  using Params = typename ExpectedValues<File>::Params;
  static constexpr auto& DefaultParams = ExpectedValues<File>::DefaultParams;
  explicit NCHG(std::shared_ptr<const File> f, const hictk::Chromosome& chrom1,
                const hictk::Chromosome& chrom2, const Params& params);
  NCHG(std::shared_ptr<const File> f, const hictk::Chromosome& chrom1,
       const hictk::Chromosome& chrom2, ExpectedValues<File> expected_values) noexcept;

  [[nodiscard]] auto params() const noexcept -> Params;

  [[nodiscard]] auto observed_matrix() const noexcept -> const ObservedMatrix<PixelIt>&;
  [[nodiscard]] auto expected_matrix() const noexcept -> const ExpectedMatrix<PixelIt>&;

  [[nodiscard]] auto compute(const hictk::GenomicInterval& range, double bad_bin_fraction) const
      -> Stats;
  [[nodiscard]] auto compute(const hictk::GenomicInterval& range1,
                             const hictk::GenomicInterval& range2, double bad_bin_fraction) const
      -> Stats;

  [[nodiscard]] auto cbegin(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> iterator;
  [[nodiscard]] auto cend(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> iterator;

  [[nodiscard]] auto begin(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> iterator;
  [[nodiscard]] auto end(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const
      -> iterator;

  [[nodiscard]] auto compute_expected_profile() const
      -> std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>;

 private:
  [[nodiscard]] static auto init_exp_matrix(const hictk::Chromosome& chrom1,
                                            const hictk::Chromosome& chrom2, const File& fp,
                                            const ExpectedValues<File>& expected_values)
      -> std::shared_ptr<const ExpectedMatrix<PixelIt>>;
  [[nodiscard]] static auto init_obs_matrix(const hictk::Chromosome& chrom1,
                                            const hictk::Chromosome& chrom2, const File& fp,
                                            const std::vector<bool>& bin1_mask,
                                            const std::vector<bool>& bin2_mask, double mad_max_,
                                            std::uint64_t min_delta_, std::uint64_t max_delta_)
      -> std::shared_ptr<const ObservedMatrix<PixelIt>>;

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

 public:
  class iterator {
    PixelIt _pixel_it{};
    PixelIt _sentinel_it{};

    std::shared_ptr<const ObservedMatrix<PixelIt>> _obs{};
    std::shared_ptr<const ExpectedMatrix<PixelIt>> _exp{};

    std::shared_ptr<const std::vector<bool>> _bin_mask1{};
    std::shared_ptr<const std::vector<bool>> _bin_mask2{};

    std::shared_ptr<std::vector<double>> _buffer{std::make_shared<std::vector<double>>()};

    std::uint64_t _min_delta{};
    std::uint64_t _max_delta{};

    mutable Stats _value{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Stats;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelIt pixel_it, PixelIt sentinel_it,
             std::shared_ptr<const ObservedMatrix<PixelIt>> obs,
             std::shared_ptr<const ExpectedMatrix<PixelIt>> exp,
             std::shared_ptr<const std::vector<bool>> bin_mask1,
             std::shared_ptr<const std::vector<bool>> bin_mask2, std::uint64_t min_delta,
             std::uint64_t max_delta) noexcept;

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
  };
};

}  // namespace nchg

#include "../../nchg_impl.hpp"
