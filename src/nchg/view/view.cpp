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

#include <fmt/compile.h>
#include <fmt/format.h>

#include <hictk/numeric_utils.hpp>

#include "nchg/config.hpp"
#include "nchg/io.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {
struct NCHGFilterResult {
  hictk::Pixel<std::uint32_t> pixel{};
  double expected{};
  double pval{};
  double pval_corrected{};
  double log_ratio{};
  double odds_ratio{};
  double omega{};
};

using NCHGComputeResult = NCHGResult;

template <typename N>
[[nodiscard]] static N parse_numeric(std::string_view tok) {
  return hictk::internal::parse_numeric_or_throw<N>(tok);
}

using ParquetFile =
    std::variant<ParquetStatsFile<NCHGComputeResult>, ParquetStatsFile<NCHGFilterResult>>;

static std::tuple<std::string, std::uint32_t, std::uint32_t> parse_ucsc(std::string buffer) {
  try {
    if (buffer.empty()) {
      throw std::runtime_error("query is empty");
    }

    const auto p1 = buffer.find_last_of(':');
    auto p2 = buffer.find_last_of('-');

    if (p1 == std::string::npos && p2 == std::string::npos) {
      return std::make_tuple(buffer, 0, std::numeric_limits<std::uint32_t>::max());
    }

    if (p1 == std::string::npos || p2 == std::string::npos || p1 > p2) {
      throw std::runtime_error(fmt::format(FMT_STRING("query \"{}\" is malformed"), buffer));
    }

    if (buffer.find(',', p1) != std::string::npos) {
      buffer.erase(std::remove(buffer.begin() + std::ptrdiff_t(p1), buffer.end(), ','),
                   buffer.end());
      p2 = buffer.find_last_of('-');
    }

    const auto chrom = buffer.substr(0, p1);
    const auto start_str = buffer.substr(p1 + 1, p2 - (p1 + 1));
    const auto end_str = buffer.substr(p2 + 1);

    return std::make_tuple(chrom, parse_numeric<std::uint32_t>(start_str),
                           parse_numeric<std::uint32_t>(end_str));
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to parse UCSC string \"{}\": {}"), buffer, e.what()));
  }
}

template <typename T, typename std::enable_if_t<std::is_same_v<T, NCHGComputeResult>>* = nullptr>
static void print_header() {
  fmt::print(
      FMT_STRING("chrom1\tstart1\tend1\t"
                 "chrom2\tstart2\tend2\t"
                 "pvalue\t"
                 "observed_count\t"
                 "expected_count\t"
                 "log_ratio\t"
                 "odds_ratio\t"
                 "omega\n"));
}

template <typename T, typename std::enable_if_t<std::is_same_v<T, NCHGFilterResult>>* = nullptr>
static void print_header() {
  fmt::print(
      FMT_STRING("chrom1\tstart1\tend1\t"
                 "chrom2\tstart2\tend2\t"
                 "pvalue\t"
                 "pvalue_corrected\t"
                 "observed_count\t"
                 "expected_count\t"
                 "log_ratio\t"
                 "odds_ratio\t"
                 "omega\n"));
}

static void process_record(const NCHGComputeResult& record, bool filter_by_coords1,
                           std::string_view chrom1, std::uint32_t start1, std::uint32_t end1,
                           bool filter_by_coords2, std::string_view chrom2, std::uint32_t start2,
                           std::uint32_t end2, double pvalue_cutoff,
                           [[maybe_unused]] double log_ratio_cutoff) {
  const auto& bin1 = record.pixel.coords.bin1;
  if (filter_by_coords1 &&
      (bin1.chrom().name() != chrom1 || bin1.start() < start1 || bin1.start() >= end1)) {
    return;
  }

  const auto& bin2 = record.pixel.coords.bin2;
  if (filter_by_coords2 &&
      (bin2.chrom().name() != chrom2 || bin2.start() < start2 || bin2.start() >= end2)) {
    return;
  }

  if (record.pval > pvalue_cutoff) {
    return;
  }

  if (std::isfinite(record.log_ratio) && record.log_ratio < log_ratio_cutoff) {
    return;
  }

  fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\t"
                         "{:s}\t{:d}\t{:d}\t"
                         "{:g}\t{:d}\t{:g}\t{:g}\t{:g}\n"),
             record.pixel.coords.bin1.chrom().name(), record.pixel.coords.bin1.start(),
             record.pixel.coords.bin1.end(), record.pixel.coords.bin2.chrom().name(),
             record.pixel.coords.bin2.start(), record.pixel.coords.bin2.end(), record.pval,
             record.pixel.count, record.expected, record.log_ratio, record.odds_ratio,
             record.omega);
}

static void process_record(const NCHGFilterResult& record, bool filter_by_coords1,
                           std::string_view chrom1, std::uint32_t start1, std::uint32_t end1,
                           bool filter_by_coords2, std::string_view chrom2, std::uint32_t start2,
                           std::uint32_t end2, double pvalue_cutoff, double log_ratio_cutoff) {
  const auto& bin1 = record.pixel.coords.bin1;
  if (filter_by_coords1 &&
      (bin1.chrom().name() != chrom1 || bin1.start() < start1 || bin1.start() >= end1)) {
    return;
  }

  const auto& bin2 = record.pixel.coords.bin2;
  if (filter_by_coords2 &&
      (bin2.chrom().name() != chrom2 || bin2.start() < start2 || bin2.start() >= end2)) {
    return;
  }

  if (record.pval_corrected > pvalue_cutoff) {
    return;
  }

  if (std::isfinite(record.log_ratio) && record.log_ratio < log_ratio_cutoff) {
    return;
  }

  fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\t"
                         "{:s}\t{:d}\t{:d}\t"
                         "{:g}\t{:g}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n"),
             record.pixel.coords.bin1.chrom().name(), record.pixel.coords.bin1.start(),
             record.pixel.coords.bin1.end(), record.pixel.coords.bin2.chrom().name(),
             record.pixel.coords.bin2.start(), record.pixel.coords.bin2.end(), record.pval,
             record.pval_corrected, record.pixel.count, record.expected, record.log_ratio,
             record.odds_ratio, record.omega);
}

int run_nchg_view(const ViewConfig& c) {
  ParquetFile f{ParquetStatsFile<NCHGResult>{}};

  if (ParquetStatsFile<NCHGResult>::is_nchg_compute_parquet(c.input_path)) {
    f = ParquetStatsFile<NCHGResult>{c.input_path};
  } else if (ParquetStatsFile<NCHGFilterResult>::is_nchg_filter_parquet(c.input_path)) {
    f = ParquetStatsFile<NCHGFilterResult>{c.input_path};
  } else {
    throw std::runtime_error(fmt::format(
        FMT_STRING("file {} does not appear to have been generated by NCHG"), c.input_path));
  }

  std::visit(
      [&](auto& ff) {
        const auto [chrom1, start1, end1] = parse_ucsc(c.range1);
        const auto [chrom2, start2, end2] = parse_ucsc(c.range2);

        const auto filter_by_coords1 = chrom1 != "all";
        const auto filter_by_coords2 = chrom2 != "all";

        if (c.with_header) {
          using T = remove_cvref_t<decltype(*ff.begin())>;
          print_header<T>();
        }

        for (const auto& record : ff) {
          process_record(record, filter_by_coords1, chrom1, start1, end1, filter_by_coords2, chrom2,
                         start2, end2, c.pvalue_cutoff, c.log_ratio_cutoff);
        }
      },
      f);

  return 0;
}
}  // namespace nchg
