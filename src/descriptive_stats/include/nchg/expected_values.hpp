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
#include <parallel_hashmap/phmap.h>

#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <highfive/H5File.hpp>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "nchg/concepts.hpp"
#include "nchg/expected_matrix.hpp"

namespace nchg {

class ExpectedValues {
  using N = double;
  std::shared_ptr<const hictk::File> _fp;
  std::uint32_t _resolution{};

  using ChromPair = std::pair<hictk::Chromosome, hictk::Chromosome>;
  using BinMask = std::shared_ptr<const std::vector<bool>>;

  std::vector<double> _expected_weights;
  phmap::btree_map<hictk::Chromosome, double> _expected_scaling_factors;
  phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>> _bin_masks;

  phmap::btree_map<ChromPair, double> _expected_values_trans;

  double _mad_max{};
  std::uint64_t _min_delta{};
  std::uint64_t _max_delta{std::numeric_limits<std::uint64_t>::max()};
  double _bin_aggregation_possible_distances_cutoff{};
  double _bin_aggregation_observed_distances_cutoff{};
  bool _interpolate{};
  double _interpolation_qtile{};
  std::uint32_t _interpolation_window_size{};
  bool _seeded{};

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
      .mad_max = 5.0,
      .min_delta = 40'000,
      .max_delta = std::numeric_limits<std::uint64_t>::max(),
      .bin_aggregation_possible_distances_cutoff = 10'000,
      .bin_aggregation_observed_distances_cutoff = 100'000,
      .interpolate = true,
      .interpolation_qtile = 0.975,
      .interpolation_window_size = 750'000
  };
  // clang-format on

  ExpectedValues() = default;
  explicit ExpectedValues(
      std::shared_ptr<const hictk::File> file, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});
  static ExpectedValues cis_only(
      std::shared_ptr<const hictk::File> file, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});
  static ExpectedValues trans_only(
      std::shared_ptr<const hictk::File> file, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});
  static ExpectedValues chromosome_pair(
      std::shared_ptr<const hictk::File> file, const hictk::Chromosome& chrom1,
      const hictk::Chromosome& chrom2, const Params& params_ = DefaultParams,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask = {});

  static ExpectedValues deserialize(const std::filesystem::path& path);

  [[nodiscard]] std::uint32_t resolution() const noexcept;
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

  [[nodiscard]] ExpectedMatrixStats expected_matrix(const hictk::Chromosome& chrom) const;
  template <typename Pixels>
    requires PixelRange<Pixels>
  [[nodiscard]] ExpectedMatrixStats expected_matrix(const hictk::Chromosome& chrom,
                                                    const hictk::BinTable& bins,
                                                    const Pixels& pixels) const;
  [[nodiscard]] ExpectedMatrixStats expected_matrix(const hictk::Chromosome& chrom1,
                                                    const hictk::Chromosome& chrom2) const;
  template <typename Pixels>
    requires PixelRange<Pixels>
  [[nodiscard]] ExpectedMatrixStats expected_matrix(const hictk::Chromosome& chrom1,
                                                    const hictk::Chromosome& chrom2,
                                                    const hictk::BinTable& bins,
                                                    const Pixels& pixels) const;

  void serialize(const std::filesystem::path& path) const;

  template <typename File>
    requires HictkSingleResFile<File>
  [[nodiscard]] static auto init_pixel_merger_cis(const File& f);

  [[nodiscard]] bool cis_only() const noexcept;
  [[nodiscard]] bool trans_only() const noexcept;

 private:
  void compute_expected_values_cis(
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask);
  void compute_expected_values_trans(
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask);

  void add_bin_mask(
      const hictk::Chromosome& chrom, std::vector<bool> mask,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask_seed);
  void add_bin_mask(
      const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
      std::pair<std::vector<bool>, std::vector<bool>> masks,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask_seed);

  static void serialize_attributes(HighFive::File& f, const Params& params);
  static void serialize_chromosomes(HighFive::File& f, const hictk::Reference& chroms);
  static void serialize_bin_masks(
      HighFive::File& f, const phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>>& bin_masks,
      bool seeded);
  static void serialize_cis_profiles(
      HighFive::File& f, const std::vector<double>& profile,
      const phmap::btree_map<hictk::Chromosome, double>& scaling_factors);
  static void serialize_trans_profiles(HighFive::File& f,
                                       const phmap::btree_map<ChromPair, double>& nnz_avg_values);

  [[nodiscard]] static auto deserialize_attributes(const HighFive::File& f) -> Params;
  [[nodiscard]] static hictk::Reference deserialize_chromosomes(const HighFive::File& f);
  [[nodiscard]] static auto deserialize_bin_masks(HighFive::File& f)
      -> phmap::btree_map<ChromPair, std::pair<BinMask, BinMask>>;
  [[nodiscard]] static std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
  deserialize_cis_profiles(const HighFive::File& f);
  [[nodiscard]] static auto deserialize_trans_profiles(const HighFive::File& f)
      -> phmap::btree_map<ChromPair, double>;

  static void merge_bin_masks(std::vector<bool>& mask1, const std::vector<bool>& mask2);

  template <typename File>
    requires HictkSingleResFile<File>
  void init_bin_masks(
      const File& f,
      const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>>& bin_mask_seed);
};
}  // namespace nchg

#include "../expected_values_impl.hpp"
