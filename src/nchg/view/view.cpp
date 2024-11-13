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

#include <concepts>
#include <hictk/genomic_interval.hpp>

#include "nchg/nchg.hpp"
#include "nchg/parquet_stats_file.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {
struct NCHGFilterResult {
  hictk::Pixel<std::uint64_t> pixel{};
  double expected{};
  double pval{};
  double pval_corrected{};
  double log_ratio{};
  double odds_ratio{};
  double omega{};
};

using NCHGComputeResult = NCHGResult;

template <typename T>
  requires std::same_as<T, NCHGComputeResult>
static void print_header() {
  fmt::print(
      "chrom1\tstart1\tend1\t"
      "chrom2\tstart2\tend2\t"
      "pvalue\t"
      "observed_count\t"
      "expected_count\t"
      "log_ratio\t"
      "odds_ratio\t"
      "omega\n");
}

template <typename T>
  requires std::same_as<T, NCHGFilterResult>
static void print_header() {
  fmt::print(
      "chrom1\tstart1\tend1\t"
      "chrom2\tstart2\tend2\t"
      "pvalue\t"
      "pvalue_corrected\t"
      "observed_count\t"
      "expected_count\t"
      "log_ratio\t"
      "odds_ratio\t"
      "omega\n");
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

template <typename T>
[[nodiscard]] static int process_records(ParquetStatsFile& f, const ViewConfig& c) {
  const auto [chrom1, start1, end1] = hictk::GenomicInterval::parse_ucsc(c.range1);
  const auto [chrom2, start2, end2] = hictk::GenomicInterval::parse_ucsc(c.range2);

  const auto filter_by_coords1 = chrom1 != "all";
  const auto filter_by_coords2 = chrom2 != "all";

  if (c.with_header) {
    print_header<T>();
  }

  auto first = f.begin<T>();
  auto last = f.end<T>();

  // apple-clang < 16 chokes if we use std::for_each()
  while (first != last) {
    process_record(*first, filter_by_coords1, chrom1, start1, end1, filter_by_coords2, chrom2,
                   start2, end2, c.pvalue_cutoff, c.log_ratio_cutoff);
    std::ignore = ++first;
  }
  return 0;
}

int run_command(const ViewConfig& c) {
  using enum ParquetStatsFile::RecordType;

  ParquetStatsFile f{c.input_path, infer};

  if (f.record_type() == NCHGCompute) {
    return process_records<NCHGComputeResult>(f, c);
  }

  assert(f.record_type() == NCHGFilter);
  return process_records<NCHGFilterResult>(f, c);
}

}  // namespace nchg
