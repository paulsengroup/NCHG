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
  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{std::numeric_limits<std::uint64_t>::max()};

  using Key = std::pair<hictk::Chromosome, hictk::Chromosome>;
  phmap::flat_hash_map<Key, std::shared_ptr<const ObservedMatrix<PixelIt>>> _obs_matrices{};
  phmap::flat_hash_map<Key, std::shared_ptr<const ExpectedMatrix<PixelIt>>> _exp_matrices{};
  ExpectedValues<File> _expected_values{};

  mutable std::vector<double> _nchg_pval_buffer{};

 public:
  explicit NCHG(std::shared_ptr<const File> f, std::uint64_t min_delta = 40'000,
                std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());
  NCHG(std::shared_ptr<const File> f, ExpectedValues<File> expected_values,
       std::uint64_t min_delta = 40'000,
       std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());

  [[nodiscard]] static NCHG cis_only(
      std::shared_ptr<const File> f, std::uint64_t min_delta = 40'000,
      std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());
  [[nodiscard]] static NCHG trans_only(std::shared_ptr<const File> f);
  [[nodiscard]] static NCHG chromosome_pair(
      std::shared_ptr<const File> f, const hictk::Chromosome& chrom1,
      const hictk::Chromosome& chrom2, std::uint64_t min_delta = 40'000,
      std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());

  [[nodiscard]] auto observed_matrix(const hictk::Chromosome& chrom) const
      -> const ObservedMatrix<PixelIt>&;
  [[nodiscard]] auto observed_matrix(const hictk::Chromosome& chrom1,
                                     const hictk::Chromosome& chrom2) const
      -> const ObservedMatrix<PixelIt>&;

  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom) const
      -> const ExpectedMatrix<PixelIt>&;
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom1,
                                     const hictk::Chromosome& chrom2) const
      -> const ExpectedMatrix<PixelIt>&;

  void init_matrices();
    void init_cis_matrices();
    void init_trans_matrices();
  void init_matrix(const hictk::Chromosome& chrom);
  void init_matrix(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2);

  void erase_matrices() noexcept;
  void erase_matrix(const hictk::Chromosome& chrom);
  void erase_matrix(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2);

  [[nodiscard]] auto compute(const hictk::GenomicInterval& range) const -> Stats;
  [[nodiscard]] auto compute(const hictk::GenomicInterval& range1,
                             const hictk::GenomicInterval& range2) const -> Stats;

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

    std::shared_ptr<const ObservedMatrix<PixelIt>> _obs{};
    std::shared_ptr<const ExpectedMatrix<PixelIt>> _exp{};

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
    iterator(PixelIt pixel_it, std::shared_ptr<const ObservedMatrix<PixelIt>> obs,
             std::shared_ptr<const ExpectedMatrix<PixelIt>> exp, std::uint64_t min_delta,
             std::uint64_t max_delta);

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
