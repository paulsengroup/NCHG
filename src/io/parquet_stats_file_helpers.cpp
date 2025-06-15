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

#include "nchg/parquet_stats_file_helpers.hpp"

#include <arrow/type.h>
#include <arrow/util/key_value_metadata.h>
#include <fmt/format.h>
#include <parquet/platform.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nchg {

std::shared_ptr<arrow::Schema> make_schema(
    std::shared_ptr<const arrow::KeyValueMetadata> metadata) {
  const auto chrom_dtype = arrow::dictionary(arrow::int32(), arrow::utf8());

  return arrow::schema(
      {
          // clang-format off
          arrow::field("chrom1",         chrom_dtype,      false),
          arrow::field("start1",         arrow::uint32(),  false),
          arrow::field("end1",           arrow::uint32(),  false),
          arrow::field("chrom2",         chrom_dtype,      false),
          arrow::field("start2",         arrow::uint32(),  false),
          arrow::field("end2",           arrow::uint32(),  false),
          arrow::field("pvalue",         arrow::float64(), false),
          arrow::field("observed_count", arrow::uint64(),  false),
          arrow::field("expected_count", arrow::float64(), false),
          arrow::field("log_ratio",      arrow::float64(), false),
          arrow::field("odds_ratio",     arrow::float64(), false),
          arrow::field("omega",          arrow::float64(), false)
          // clang-format on
      },
      std::move(metadata));
}

std::shared_ptr<arrow::Schema> make_schema_with_padj(
    std::shared_ptr<const arrow::KeyValueMetadata> metadata) {
  const auto chrom_dtype = arrow::dictionary(arrow::int32(), arrow::utf8());

  return arrow::schema(
      {
          // clang-format off
          arrow::field("chrom1",           chrom_dtype,      false),
          arrow::field("start1",           arrow::uint32(),  false),
          arrow::field("end1",             arrow::uint32(),  false),
          arrow::field("chrom2",           chrom_dtype,      false),
          arrow::field("start2",           arrow::uint32(),  false),
          arrow::field("end2",             arrow::uint32(),  false),
          arrow::field("pvalue",           arrow::float64(), false),
          arrow::field("pvalue_corrected", arrow::float64(), false),
          arrow::field("observed_count",   arrow::uint64(),  false),
          arrow::field("expected_count",   arrow::float64(), false),
          arrow::field("log_ratio",        arrow::float64(), false),
          arrow::field("odds_ratio",       arrow::float64(), false),
          arrow::field("omega",            arrow::float64(), false)
          // clang-format on
      },
      std::move(metadata));
}

parquet::Compression::type parse_parquet_compression(std::string_view method) {
  if (method == "zstd") {
    return parquet::Compression::ZSTD;
  }
  if (method == "lz4") {
    return parquet::Compression::LZ4;
  }
  throw std::runtime_error(fmt::format("unrecognized compression method \"{}\"", method));
}

}  // namespace nchg
