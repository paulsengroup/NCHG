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

#include <cassert>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <memory>
#include <type_traits>
#include <variant>

#include "nchg/common.hpp"
#include "nchg/config.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/tools.hpp"

namespace nchg {

template <typename FilePtr>
static void process_all_chromosomes(FilePtr f, const ExpectedConfig &c) {
  const ExpectedValues evs(
      f, {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
          c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
          c.interpolation_qtile, c.interpolation_window_size});
  if (c.force) {
    std::filesystem::remove(c.output_path);
  }
  evs.serialize(c.output_path);
}

template <typename FilePtr>
static void process_cis_chromosomes(FilePtr f, const ExpectedConfig &c) {
  using File = remove_cvref_t<decltype(*f)>;
  const auto evs = ExpectedValues<File>::cis_only(
      f, {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
          c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
          c.interpolation_qtile, c.interpolation_window_size});
  if (c.force) {
    std::filesystem::remove(c.output_path);
  }
  evs.serialize(c.output_path);
}

template <typename FilePtr>
static void process_trans_chromosomes(FilePtr f, const ExpectedConfig &c) {
  using File = remove_cvref_t<decltype(*f)>;
  const auto evs = ExpectedValues<File>::trans_only(
      f, {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
          c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
          c.interpolation_qtile, c.interpolation_window_size});
  if (c.force) {
    std::filesystem::remove(c.output_path);
  }
  evs.serialize(c.output_path);
}

template <typename FilePtr>
static void process_one_chromosome_pair(FilePtr f, const ExpectedConfig &c) {
  assert(c.chrom1 != "all");
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  using File = remove_cvref_t<decltype(*std::declval<FilePtr>())>;
  const auto evs = ExpectedValues<File>::chromosome_pair(
      f, chrom1, chrom2,
      {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
       c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
       c.interpolation_qtile, c.interpolation_window_size});
  if (c.force) {
    std::filesystem::remove(c.output_path);
  }
  evs.serialize(c.output_path);
}

// clang-format off
  using HiCFilePtr =
      std::variant<
          std::shared_ptr<const hictk::cooler::File>,
          std::shared_ptr<const hictk::hic::File>>;
// clang-format on

[[nodiscard]] static HiCFilePtr open_file_ptr(const std::filesystem::path &path,
                                              std::uint32_t resolution) {
  hictk::File f_(path.string(), resolution);
  return {std::visit(
      [&](auto &&ff) {
        using FileT = std::remove_reference_t<decltype(ff)>;
        return HiCFilePtr{std::make_shared<const FileT>(std::forward<FileT>(ff))};
      },
      f_.get())};
}

int run_nchg_expected(const ExpectedConfig &c) {
  const auto f = open_file_ptr(c.input_path, c.resolution);

  std::visit(
      [&](const auto &f_) {
        if (c.cis_only) {
          process_cis_chromosomes(f_, c);
          return;
        }
        if (c.trans_only) {
          process_trans_chromosomes(f_, c);
          return;
        }
        if (c.chrom1 == "all") {
          assert(c.chrom2 == "all");
          process_all_chromosomes(f_, c);
          return;
        }
        process_one_chromosome_pair(f_, c);
      },
      f);

  return 0;
}

}  // namespace nchg
