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

#include "nchg/parquet_helpers.hpp"

#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/type.h>
#include <arrow/util/key_value_metadata.h>
#include <parquet/platform.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>

#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nchg {

static std::shared_ptr<arrow::KeyValueMetadata> generate_schema_metadata(
    const hictk::Reference &chroms) {
  std::vector<std::string> keys{};
  std::vector<std::string> values{};
  for (const auto &chrom : chroms) {
    if (!chrom.is_all()) {
      keys.emplace_back(chrom.name());
      values.emplace_back(fmt::to_string(chrom.size()));
    }
  }

  return std::make_shared<arrow::KeyValueMetadata>(std::move(keys), std::move(values));
}

std::shared_ptr<arrow::Schema> get_schema(const hictk::Reference &chroms) {
  const auto chrom_dtype = arrow::dictionary(arrow::int32(), arrow::utf8());
  const auto metadata = generate_schema_metadata(chroms);

  return arrow::schema({
      // clang-format off
      arrow::field("chrom1",         chrom_dtype,          false, metadata),
      arrow::field("start1",         arrow::uint32(),  false),
      arrow::field("end1",           arrow::uint32(),  false),
      arrow::field("chrom2",         chrom_dtype,          false, metadata),
      arrow::field("start2",         arrow::uint32(),  false),
      arrow::field("end2",           arrow::uint32(),  false),
      arrow::field("pvalue",         arrow::float64(), false),
      arrow::field("observed_count", arrow::uint64(),  false),
      arrow::field("expected_count", arrow::float64(), false),
      arrow::field("log_ratio",      arrow::float64(), false),
      arrow::field("odds_ratio",     arrow::float64(), false),
      arrow::field("omega",          arrow::float64(), false)
      // clang-format on
  });
}

std::shared_ptr<arrow::Schema> get_schema_padj(const hictk::Reference &chroms) {
  const auto chrom_dtype = arrow::dictionary(arrow::int32(), arrow::utf8());
  const auto metadata = generate_schema_metadata(chroms);

  return arrow::schema({
      // clang-format off
      arrow::field("chrom1",           chrom_dtype,          false, metadata),
      arrow::field("start1",           arrow::uint32(),  false),
      arrow::field("end1",             arrow::uint32(),  false),
      arrow::field("chrom2",           chrom_dtype,          false, metadata),
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
  });
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
