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
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <variant>

namespace nchg {

// NOLINTBEGIN(clang-analyzer-optin.performance.Padding)
struct CartesianProductConfig {
  std::filesystem::path path_to_domains;
  std::filesystem::path path_to_chrom_sizes;

  std::optional<std::string> chrom1;
  std::optional<std::string> chrom2;

  bool process_cis{true};
  bool process_trans{true};

  std::uint8_t verbosity{3};
};

struct ChecksumConfig {
  std::vector<std::filesystem::path> files;

  std::uint8_t verbosity{3};
};

struct ComputePvalConfig {
  enum class DomainAggregationStrategy : std::uint8_t { AUTO, SINGLE_PASS, MULTI_PASS };

  std::filesystem::path exec;
  std::filesystem::path path_to_hic;
  std::filesystem::path output_prefix;
  std::filesystem::path output_path;
  bool force{false};

  std::optional<std::uint32_t> resolution;
  std::optional<std::string> chrom1;
  std::optional<std::string> chrom2;
  std::filesystem::path path_to_expected_values;
  std::filesystem::path path_to_domains;

  bool skip_empty_matrices{true};
  bool compute_cis{true};
  bool compute_trans{true};
  DomainAggregationStrategy domain_aggregation_stategy{DomainAggregationStrategy::AUTO};

  double mad_max{5.0};
  std::uint64_t min_delta{40'000};
  std::uint64_t max_delta{std::numeric_limits<std::uint64_t>::max()};
  double bin_aggregation_possible_distances_cutoff{0};
  double bin_aggregation_observed_distances_cutoff{5'000};
  bool interpolate_expected_values{true};
  double interpolation_qtile{0.975};
  std::uint32_t interpolation_window_size{750'000};
  double bad_bin_fraction{0.1};
  std::filesystem::path path_to_bin_mask;

  std::string compression_method{"zstd"};
  std::uint8_t compression_lvl{9};
  std::size_t threads{1};

  std::uint8_t verbosity{3};
  std::string log_message_queue;
};

struct FilterConfig {
  std::filesystem::path input_path;
  std::filesystem::path output_path;
  bool force{false};

  double fdr{0.01};
  double log_ratio{2.0};
  bool drop_non_significant{true};

  bool correct_cis_trans_separately{true};
  bool correct_chrom_chrom_separately{false};

  std::size_t threads{2};
  std::string compression_method{"zstd"};
  std::uint8_t compression_lvl{9};

  std::uint8_t verbosity{3};
};

struct ExpectedConfig {
  std::filesystem::path input_path;
  std::optional<std::uint32_t> resolution;

  std::filesystem::path output_path;
  bool force{false};

  std::optional<std::string> chrom1;
  std::optional<std::string> chrom2;
  bool cis_only{false};
  bool trans_only{false};

  double mad_max{5.0};
  std::uint64_t min_delta{40'000};
  std::uint64_t max_delta{std::numeric_limits<std::uint64_t>::max()};
  double bin_aggregation_possible_distances_cutoff{1'000};
  double bin_aggregation_observed_distances_cutoff{100'000};
  bool interpolate_expected_values{true};
  double interpolation_qtile{0.975};
  std::uint32_t interpolation_window_size{750'000};
  std::filesystem::path path_to_bin_mask;

  std::uint8_t verbosity{3};
};

struct MergeConfig {
  std::filesystem::path input_prefix;
  std::filesystem::path output_path;
  bool force{false};
  bool ignore_report_file{false};

  std::size_t threads{2};
  std::string compression_method{"zstd"};
  std::uint8_t compression_lvl{9};

  std::uint8_t verbosity{3};
};

struct ViewConfig {
  std::filesystem::path input_path;

  std::string range1{"all"};
  std::string range2{"all"};

  double pvalue_cutoff{1.0};
  double log_ratio_cutoff{-std::numeric_limits<double>::infinity()};

  bool with_header{true};

  std::uint8_t verbosity{3};
};

struct MetadataConfig {
  std::filesystem::path input_path;

  bool raw{false};

  std::uint8_t verbosity{0};
};

// NOLINTEND(clang-analyzer-optin.performance.Padding)

// clang-format off
using Config = std::variant<
    std::monostate,
    CartesianProductConfig,
    ChecksumConfig,
    ComputePvalConfig,
    ExpectedConfig,
    FilterConfig,
    MergeConfig,
    MetadataConfig,
    ViewConfig>;
// clang-format on

}  // namespace nchg
