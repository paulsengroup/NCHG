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

#include "nchg/concepts.hpp"
#include "nchg/expected_matrix.hpp"

namespace nchg {

template <typename File>
  requires HictkSingleResFile<File>
class ExpectedValues {
  using N = std::uint32_t;

  std::shared_ptr<const File> _fp{};
  std::uint32_t _resolution{};

  using ChromPair = std::pair<hictk::Chromosome, hictk::Chromosome>;
  using BinMask = std::shared_ptr<const std::vector<bool>>;

  std::vector<double> _expected_weights{};
  phmap::btree_map<hictk::Chromosome, double> _expected_scaling_factors{};
  phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>> _bin_masks{};

  phmap::btree_map<ChromPair, double> _expected_values_trans{};

  double _mad_max{};
  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{std::numeric_limits<std::uint64_t>::max()};
  double _bin_aggregation_possible_distances_cutoff{};
  double _bin_aggregation_observed_distances_cutoff{};
  bool _interpolate{};
  double _interpolation_qtile{};
  std::uint32_t _interpolation_window_size{};

 public:
  struct Params {
    double mad_max{};
    std::uint64_t min_delta{};
    std::uint64_t max_delta{};
    double bin_aggregation_possible_distances_cutoff{};
    double bin_aggregation_observed_distances_cutoff{};
    bool interpolate{};
    double interpolation_qtile{};
    std::uint32_t interpolation_window_size{};
  };

  // clang-format off
  static constexpr Params DefaultParams{
      5.0,
      40'000,
      std::numeric_limits<std::uint64_t>::max(),
      10'000,
      100'000,
      true,
      0.975,
      750'000
  };
  // clang-format on

  template <typename OtherFile>
    requires HictkSingleResFile<OtherFile>
  friend class ExpectedValues;

  explicit ExpectedValues(
      std::shared_ptr<const File> file, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});
  static ExpectedValues cis_only(
      std::shared_ptr<const File> file, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});
  static ExpectedValues trans_only(
      std::shared_ptr<const File> file, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});
  static ExpectedValues chromosome_pair(
      std::shared_ptr<const File> file, const hictk::Chromosome& chrom1,
      const hictk::Chromosome& chrom2, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});

  static ExpectedValues deserialize(const std::filesystem::path& path);

  template <typename OutFile>
    requires HictkSingleResFile<OutFile>
  [[nodiscard]] ExpectedValues<OutFile> cast() const;

  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] const std::vector<double>& weights() const noexcept;
  [[nodiscard]] auto params() const noexcept -> Params;

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

  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom) const;
  template <typename Pixels>
    requires PixelRange<Pixels>
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom, const hictk::BinTable& bins,
                                     const Pixels& pixels) const;
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom1,
                                     const hictk::Chromosome& chrom2) const;
  template <typename Pixels>
    requires PixelRange<Pixels>
  [[nodiscard]] auto expected_matrix(const hictk::Chromosome& chrom1,
                                     const hictk::Chromosome& chrom2, const hictk::BinTable& bins,
                                     const Pixels& pixels) const;

  void serialize(const std::filesystem::path& path) const;

 private:
  void compute_expected_values_cis(
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask);
  void compute_expected_values_trans(
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask);

  void add_bin_mask(
      const hictk::Chromosome& chrom, std::vector<bool>&& mask,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask_seed);
  void add_bin_mask(
      const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
      std::pair<std::vector<bool>, std::vector<bool>>&& masks2,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask_seed);

  static void serialize_attributes(HighFive::File& f, const Params& params);
  static void serialize_chromosomes(HighFive::File& f, const hictk::Reference& chroms);
  static void serialize_bin_masks(
      HighFive::File& f, const phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>>& bin_masks);
  static void serialize_cis_profiles(
      HighFive::File& f, const std::vector<double>& profile,
      const phmap::btree_map<hictk::Chromosome, double>& scaling_factors);
  static void serialize_trans_profiles(HighFive::File& f,
                                       const phmap::btree_map<ChromPair, double>& nnz_avg_values);

  [[nodiscard]] static auto deserialize_attributes(const HighFive::File& f) -> Params;
  [[nodiscard]] static hictk::Reference deserialize_chromosomes(const HighFive::File& f);
  [[nodiscard]] static auto deserialize_bin_masks(HighFive::File& f)
      -> const phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>>;
  [[nodiscard]] static std::pair<const std::vector<double>,
                                 const phmap::btree_map<hictk::Chromosome, double>>
  deserialize_cis_profiles(const HighFive::File& f);
  [[nodiscard]] static auto deserialize_trans_profiles(const HighFive::File& f)
      -> phmap::btree_map<ChromPair, double>;

  static void merge_bin_masks(std::vector<bool>& mask1, const std::vector<bool>& mask2);
};
}  // namespace nchg

#include "../expected_values_impl.hpp"
