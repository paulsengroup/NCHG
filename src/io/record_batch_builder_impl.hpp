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

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/array/array_base.h>
#include <arrow/io/file.h>
#include <arrow/util/thread_pool.h>
#include <parquet/properties.h>
#include <parquet/arrow/writer.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include "nchg/parquet_helpers.hpp"

namespace nchg {

template <typename Stats>
inline void RecordBatchBuilder::append(const Stats &s) {
  append(_chrom1, std::string{s.pixel.coords.bin1.chrom().name()});
  append(_start1, s.pixel.coords.bin1.start());
  append(_end1, s.pixel.coords.bin1.end());

  append(_chrom2, std::string{s.pixel.coords.bin2.chrom().name()});
  append(_start2, s.pixel.coords.bin2.start());
  append(_end2, s.pixel.coords.bin2.end());

  append(_pvalue, s.pval);
  append(_observed, s.pixel.count);
  append(_expected, s.expected);
  append(_log_ratio, s.log_ratio);
  append(_odds, s.odds_ratio);
  append(_omega, s.omega);

  if constexpr (has_pval_corrected<Stats>::value) {
    append(_pvalue_corrected, s.pval_corrected);
  }

  ++_i;
}

template <typename ArrayBuilder, typename T>
inline void RecordBatchBuilder::append(ArrayBuilder &builder, const T &data) {
  const auto status = builder.Append(data);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

template <typename ArrayBuilder>
inline std::shared_ptr<arrow::Array> RecordBatchBuilder::finish(ArrayBuilder &builder) {
  auto result = builder.Finish();
  if (!result.status().ok()) {
    throw std::runtime_error(result.status().ToString());
  }

  return result.MoveValueUnsafe();
}

template <typename Record>
inline std::unique_ptr<parquet::arrow::FileWriter> init_parquet_file_writer(
    const hictk::Reference &chroms, const std::filesystem::path &path, bool force,
    std::string_view compression_method, std::uint8_t compression_lvl, std::size_t threads) {
  if (path.empty()) {
    return {};
  }

  const auto schema =
      has_pval_corrected<Record>::value ? *get_schema_padj(chroms) : *get_schema(chroms);

  auto builder = parquet::WriterProperties::Builder()
                     .created_by("NCHG v0.0.2")
                     ->version(parquet::ParquetVersion::PARQUET_2_6)
                     ->data_page_version(parquet::ParquetDataPageVersion::V2)
                     ->compression(parse_parquet_compression(compression_method))
                     ->compression_level(compression_lvl)
                     ->disable_statistics()
                     ->build();

  if (threads > 1) {
    auto status = arrow::SetCpuThreadPoolCapacity(static_cast<std::int32_t>(threads - 1));
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
    status = arrow::io::SetIOThreadPoolCapacity(static_cast<std::int32_t>(threads - 1));
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  }

  auto arrow_properties = parquet::ArrowWriterProperties::Builder()
                              .set_use_threads(threads > 1)
                              ->store_schema()
                              ->build();

  const auto output_dir = path.parent_path();
  if (!output_dir.empty() && !std::filesystem::exists(output_dir)) {
    std::filesystem::create_directories(output_dir);
  }

  if (force) {
    std::filesystem::remove(path);  // NOLINT
  } else if (std::filesystem::exists(path)) {
    throw std::runtime_error(
        fmt::format("refusing to overwrite output file \"{}\". Pass --force to overwrite.", path));
  }

  SPDLOG_DEBUG("initializing file \"{}\"...", path);

  std::shared_ptr<arrow::io::FileOutputStream> f{};
  PARQUET_ASSIGN_OR_THROW(f, arrow::io::FileOutputStream::Open(path));

  auto result = parquet::arrow::FileWriter::Open(schema, arrow::default_memory_pool(), f, builder,
                                                 arrow_properties);
  if (!result.status().ok()) {
    throw std::runtime_error(result.status().ToString());
  }

  return result.MoveValueUnsafe();
}

}  // namespace nchg
