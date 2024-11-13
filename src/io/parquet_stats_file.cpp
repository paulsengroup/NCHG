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

#include "nchg/parquet_stats_file.hpp"

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/type.h>
#include <arrow/io/file.h>
#include <arrow/util/key_value_metadata.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>
#include <parquet/stream_reader.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <hictk/numeric_utils.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "nchg/common.hpp"

namespace nchg {

[[nodiscard]] static std::shared_ptr<arrow::io::ReadableFile> open_parquet_file(
    const std::filesystem::path &path) {
  try {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path));
    return fp;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to open file \"{}\" for reading: {}", path, e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to open file \"{}\" for reading: unknown error", path));
  }
}

[[nodiscard]] std::shared_ptr<arrow::Schema> get_file_schema(
    const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  try {
    auto props = parquet::default_reader_properties();
    props.set_buffer_size(1);

    std::unique_ptr<parquet::arrow::FileReader> reader;
    auto status = parquet::arrow::OpenFile(fp, arrow::default_memory_pool(), &reader);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    std::unique_ptr<arrow::RecordBatchReader> batch_reader{};
    status = reader->GetRecordBatchReader(&batch_reader);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    return batch_reader->schema();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to read file schema: {}", e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read file schema: unknown error");
  }
}

[[nodiscard]] static std::shared_ptr<const hictk::Reference> import_chromosomes_from_parquet(
    const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  try {
    auto schema = get_file_schema(fp);

    const auto metadata = schema->field(0)->metadata();
    const auto chrom_names = metadata->keys();
    std::vector<std::uint32_t> chrom_sizes(chrom_names.size(), 0);
    std::ranges::transform(metadata->values(), chrom_sizes.begin(), [](const auto &size) {
      return hictk::internal::parse_numeric_or_throw<std::uint32_t>(size);
    });

    return std::make_shared<const hictk::Reference>(chrom_names.begin(), chrom_names.end(),
                                                    chrom_sizes.begin());

  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to read chromosomes: {}", e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read chromosomes: unknown error");
  }
}

[[nodiscard]] static std::shared_ptr<parquet::StreamReader> init_parquet_stream_reader(
    std::shared_ptr<arrow::io::ReadableFile> fp, std::size_t buffer_size) {
  assert(buffer_size != 0);
  auto props = parquet::default_reader_properties();
  props.set_buffer_size(static_cast<std::int64_t>(buffer_size));

  return std::make_shared<parquet::StreamReader>(
      parquet::ParquetFileReader::Open(std::move(fp), props));
}

[[nodiscard]] static auto validate_record_type(const std::filesystem::path &path,
                                               const std::shared_ptr<arrow::io::ReadableFile> &fp,
                                               ParquetStatsFile::RecordType expected_type)
    -> ParquetStatsFile::RecordType {
  static constexpr auto expected_columns_nchg_compute = std::to_array<std::string_view>({
      "chrom1",
      "start1",
      "end1",
      "chrom2",
      "start2",
      "end2",
      "pvalue",
      "observed_count",
      "expected_count",
      "log_ratio",
      "odds_ratio",
      "omega",
  });

  static constexpr auto expected_columns_nchg_filter = std::to_array<std::string_view>({
      "chrom1",
      "start1",
      "end1",
      "chrom2",
      "start2",
      "end2",
      "pvalue",
      "pvalue_corrected",
      "observed_count",
      "expected_count",
      "log_ratio",
      "odds_ratio",
      "omega",
  });

  try {
    const auto col_names = get_file_schema(fp)->field_names();

    auto file_has_nchg_compute_records = [&] {
      return std::ranges::equal(expected_columns_nchg_compute, col_names);
    };

    auto file_has_nchg_filter_records = [&] {
      return std::ranges::equal(expected_columns_nchg_filter, col_names);
    };

    switch (expected_type) {
      using enum ParquetStatsFile::RecordType;
      case NCHGCompute: {
        if (file_has_nchg_compute_records()) {
          return NCHGCompute;
        }
        if (file_has_nchg_filter_records()) {
          throw std::runtime_error(
              "unexpected record type: expected NCHGCompute records, found NCHGFilter");
        }
        throw std::runtime_error("unexpected record type: file was likely not generated by NCHG");
      }
      case NCHGFilter: {
        if (file_has_nchg_filter_records()) {
          return NCHGFilter;
        }
        if (file_has_nchg_compute_records()) {
          throw std::runtime_error(
              "unexpected record type: expected NCHGFilter records, found NCHGCompute");
        }
        throw std::runtime_error("unexpected record type: file was likely not generated by NCHG");
      }
      case infer: {
        if (file_has_nchg_compute_records()) {
          return NCHGCompute;
        }
        if (file_has_nchg_filter_records()) {
          return NCHGFilter;
        }
        throw std::runtime_error("unexpected record type: file was likely not generated by NCHG");
      }
      default:
        unreachable_code();
    }

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("unable to detect type of file \"{}\": {}", path, e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("unable to detect type of file \"{}\": unknown error", path));
  }
}

ParquetStatsFile::ParquetStatsFile(const std::filesystem::path &path,
                                   const std::shared_ptr<arrow::io::ReadableFile> &fp,
                                   RecordType record_type, std::size_t buffer_size)
    : ParquetStatsFile(path, fp, import_chromosomes_from_parquet(fp), record_type, buffer_size) {}

ParquetStatsFile::ParquetStatsFile(const std::filesystem::path &path,
                                   std::shared_ptr<arrow::io::ReadableFile> fp,
                                   std::shared_ptr<const hictk::Reference> chromosomes,
                                   RecordType record_type, std::size_t buffer_size)
    : _type(validate_record_type(path, fp, record_type)),
      _chroms(std::move(chromosomes)),
      _sr(init_parquet_stream_reader(std::move(fp), buffer_size)) {}

ParquetStatsFile::ParquetStatsFile(const std::filesystem::path &path, RecordType record_type,
                                   std::size_t buffer_size)
    : ParquetStatsFile(path, open_parquet_file(path), record_type, buffer_size) {}

auto ParquetStatsFile::record_type() const noexcept -> RecordType { return _type; }

std::shared_ptr<const hictk::Reference> ParquetStatsFile::chromosomes() const noexcept {
  return _chroms;
}

}  // namespace nchg
