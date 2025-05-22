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

#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <concepts>
#include <filesystem>
#include <hictk/file.hpp>
#include <memory>
#include <type_traits>

#include "nchg/common.hpp"
#include "nchg/concepts.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

template <typename FPtr>
concept HictkSingleResFilePtr = requires(FPtr fp) {
  requires SmartPtr<FPtr>;
  requires HictkSingleResFile<std::remove_cvref_t<decltype(*fp)>>;
};

template <typename FilePtr>
  requires HictkSingleResFilePtr<FilePtr>
static void process_all_chromosomes(const FilePtr &f, const ExpectedConfig &c) {
  const auto mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);
  const ExpectedValues evs(
      f,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      mask);
  if (c.force) {
    std::filesystem::remove(c.output_path);  // NOLINT
  }
  evs.serialize(c.output_path);
}

template <typename FilePtr>
  requires HictkSingleResFilePtr<FilePtr>
static void process_cis_chromosomes(const FilePtr &f, const ExpectedConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();
  const auto mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  const auto evs = ExpectedValues::cis_only(
      f,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      mask);
  if (c.force) {
    std::filesystem::remove(c.output_path);  // NOLINT
  }
  evs.serialize(c.output_path);
  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("expected values have been written to \"{}\". Computation took {}", c.output_path,
              format_duration(t1 - t0));
}

template <typename FilePtr>
  requires HictkSingleResFilePtr<FilePtr>
static void process_trans_chromosomes(const FilePtr &f, const ExpectedConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();
  const auto mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  const auto evs = ExpectedValues::trans_only(
      f,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      mask);
  if (c.force) {
    std::filesystem::remove(c.output_path);  // NOLINT
  }
  evs.serialize(c.output_path);
  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("expected values have been written to \"{}\". Computation took {}", c.output_path,
              format_duration(t1 - t0));
}

template <typename FilePtr>
  requires HictkSingleResFilePtr<FilePtr>
static void process_one_chromosome_pair(const FilePtr &f, const ExpectedConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();
  const auto mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  assert(c.chrom1 != "all");
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  const auto evs = ExpectedValues::chromosome_pair(
      f, chrom1, chrom2,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      mask);
  if (c.force) {
    std::filesystem::remove(c.output_path);  // NOLINT
  }
  evs.serialize(c.output_path);
  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("expected values have been written to \"{}\". Computation took {}", c.output_path,
              format_duration(t1 - t0));
}

int run_command(const ExpectedConfig &c) {
  const auto f = std::make_shared<const hictk::File>(c.input_path.string(), c.resolution);

  if (c.cis_only) {
    process_cis_chromosomes(f, c);
    return 0;
  }
  if (c.trans_only) {
    process_trans_chromosomes(f, c);
    return 0;
  }
  if (c.chrom1 == "all") {
    assert(c.chrom2 == "all");
    process_all_chromosomes(f, c);
    return 0;
  }
  process_one_chromosome_pair(f, c);

  return 0;
}

}  // namespace nchg
