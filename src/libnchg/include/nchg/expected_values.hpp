// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>

#pragma once

#include <parallel_hashmap/btree.h>

#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <highfive/H5File.hpp>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "nchg/expected_matrix.hpp"

namespace nchg {

template <typename File>
class ExpectedValues {
  using N = std::uint32_t;
  using ThinPixelIt = decltype(std::declval<File>().fetch("chr1").template begin<N>());
  using PixelIt =
      decltype(std::declval<hictk::transformers::JoinGenomicCoords<ThinPixelIt>>().begin());
  std::shared_ptr<const File> _fp{};

  using ChromPair = std::pair<hictk::Chromosome, hictk::Chromosome>;
  using BinMask = std::shared_ptr<const std::vector<bool>>;

  std::vector<double> _expected_weights{};
  phmap::btree_map<hictk::Chromosome, double> _expected_scaling_factors{};
  phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>> _bin_masks{};

  phmap::btree_map<ChromPair, double> _expected_values_trans{};

  double _mad_max{};
  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{std::numeric_limits<std::uint64_t>::max()};

 public:
  explicit ExpectedValues(std::shared_ptr<const File> file, double mad_max = 0.0,
                          std::uint64_t min_delta = 40'000,
                          std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());
  static ExpectedValues cis_only(
      std::shared_ptr<const File> file, double mad_max = 0.0, std::uint64_t min_delta = 40'000,
      std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());
  static ExpectedValues trans_only(std::shared_ptr<const File> file, double mad_max);
  static ExpectedValues chromosome_pair(
      std::shared_ptr<const File> file, const hictk::Chromosome& chrom1,
      const hictk::Chromosome& chrom2, double mad_max = 0.0, std::uint64_t min_delta = 40'000,
      std::uint64_t max_delta = std::numeric_limits<std::uint64_t>::max());

  static ExpectedValues deserialize(const std::filesystem::path& path);

  [[nodiscard]] const std::vector<double>& weights() const noexcept;
  [[nodiscard]] double mad_max() const noexcept;
  [[nodiscard]] std::uint64_t min_delta() const noexcept;
  [[nodiscard]] std::uint64_t max_delta() const noexcept;

  [[nodiscard]] std::shared_ptr<const std::vector<bool>> bin_mask(
      const hictk::Chromosome& chrom) const;
  [[nodiscard]] std::pair<std::shared_ptr<const std::vector<bool>>,
                          std::shared_ptr<const std::vector<bool>>>
  bin_mask(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const;

  [[nodiscard]] std::vector<double> expected_values(const hictk::Chromosome& chrom,
                                                    bool rescale = true) const;
  [[nodiscard]] double expected_value(const hictk::Chromosome& chrom1,
                                      const hictk::Chromosome& chrom2) const;

  [[nodiscard]] const phmap::btree_map<hictk::Chromosome, double>& scaling_factors() const noexcept;
  [[nodiscard]] double scaling_factor(const hictk::Chromosome& chrom) const;

  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom) const
      -> ExpectedMatrix<PixelIt>;
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom, const hictk::BinTable& bins,
                                     PixelIt first_pixel, PixelIt last_pixel) const
      -> ExpectedMatrix<PixelIt>;
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom1,
                                     const hictk::Chromosome& chrom2) const
      -> ExpectedMatrix<PixelIt>;
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom1,
                                     const hictk::Chromosome& chrom2, const hictk::BinTable& bins,
                                     PixelIt first_pixel, PixelIt last_pixel) const
      -> ExpectedMatrix<PixelIt>;

  void serialize(const std::filesystem::path& path) const;

 private:
  void compute_expected_values_cis();
  void compute_expected_values_trans();

  void add_bin_mask(const hictk::Chromosome& chrom, std::vector<bool>&& mask);
  void add_bin_mask(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
                    std::pair<std::vector<bool>, std::vector<bool>>&& masks2);

  static void serialize_attributes(HighFive::File& f, double mad_max_, std::uint64_t min_delta_,
                                   std::uint64_t max_delta_);
  static void serialize_chromosomes(HighFive::File& f, const hictk::Reference& chroms);
  static void serialize_bin_masks(
      HighFive::File& f, const phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>>& bin_masks);
  static void serialize_cis_profiles(
      HighFive::File& f, const std::vector<double>& profile,
      const phmap::btree_map<hictk::Chromosome, double>& scaling_factors);
  static void serialize_trans_profiles(HighFive::File& f,
                                       const phmap::btree_map<ChromPair, double>& nnz_avg_values);

  [[nodiscard]] static std::tuple<double, std::uint64_t, std::uint64_t> deserialize_attributes(
      const HighFive::File& f);
  [[nodiscard]] static hictk::Reference deserialize_chromosomes(const HighFive::File& f);
  [[nodiscard]] static auto deserialize_bin_masks(HighFive::File& f)
      -> const phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>>;
  [[nodiscard]] static std::pair<const std::vector<double>,
                                 const phmap::btree_map<hictk::Chromosome, double>>
  deserialize_cis_profiles(const HighFive::File& f);
  [[nodiscard]] static auto deserialize_trans_profiles(const HighFive::File& f)
      -> phmap::btree_map<ChromPair, double>;
};
}  // namespace nchg

#include "../expected_values_impl.hpp"
