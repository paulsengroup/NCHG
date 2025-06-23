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

#include "nchg/parquet_stats_file_reader.hpp"

#include <arrow/io/file.h>
#include <arrow/record_batch.h>
#include <arrow/type.h>
#include <arrow/util/base64.h>
#include <arrow/util/key_value_metadata.h>
#include <fmt/format.h>
#include <parquet/arrow/reader.h>
#include <parquet/file_reader.h>
#include <parquet/stream_reader.h>
#include <zstd.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <glaze/json/json_t.hpp>
#include <hictk/numeric_utils.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/metadata.hpp"

namespace nchg {

template <typename N>
[[nodiscard]] static N parse_numeric(std::string_view tok) {
  return hictk::internal::parse_numeric_or_throw<N>(tok);
}

[[nodiscard]] static std::shared_ptr<arrow::io::ReadableFile> open_parquet_file(
    const std::filesystem::path &path) {
  try {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path))
    return fp;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to open file \"{}\" for reading: {}", path, e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to open file \"{}\" for reading: unknown error", path));
  }
}

[[nodiscard]] static std::unique_ptr<parquet::arrow::FileReader> open_arrow_file(
    const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  auto result = parquet::arrow::OpenFile(fp, arrow::default_memory_pool());
  if (!result.ok()) {
    throw std::runtime_error(result.status().ToString());
  }

  return result.MoveValueUnsafe();
}

[[nodiscard]] static std::shared_ptr<arrow::Schema> get_file_schema(
    const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  try {
    auto props = parquet::default_reader_properties();
    props.set_buffer_size(1);

    auto reader = open_arrow_file(fp);

    auto result = reader->GetRecordBatchReader();
    if (!result.ok()) {
      throw std::runtime_error(result.status().ToString());
    }

    return result.ValueUnsafe()->schema();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to read file schema: {}", e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read file schema: unknown error");
  }
}

[[nodiscard]] static glz::json_t strip_keys_from_json(glz::json_t json,
                                                      const std::vector<std::string> &keys) {
  auto &map = json.get<glz::json_t::object_t>();
  for (const auto &k : keys) {
    if (map.contains(k)) {
      map.erase(k);
    }
  }

  return json;
}

[[nodiscard]] static std::string read_attribute_or_throw(
    const std::shared_ptr<const arrow::Schema> &schema, std::string_view key) {
  const auto &metadata = schema->metadata();
  int i = -1;

  if (metadata) {
    i = schema->metadata()->FindKey(key);
  }

  if (i == -1) {
    throw std::out_of_range(fmt::format("failed to read {} attribute: key not found", key));
  }

  return schema->metadata()->value(i);
}

[[nodiscard]] static std::string read_metadata_or_throw(
    const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  try {
    const auto schema = get_file_schema(fp);

    const auto compression = read_attribute_or_throw(schema, "NCHG:metadata-compression");
    auto metadata = read_attribute_or_throw(schema, "NCHG:metadata");

    if (compression == "None") {
      return metadata;
    }

    if (compression == "zstd") {
      const auto buffer_size =
          parse_numeric<std::size_t>(read_attribute_or_throw(schema, "NCHG:metadata-size"));
      std::string buffer(buffer_size, '\0');
      metadata = arrow::util::base64_decode(std::string{metadata});
      const auto decompressed_size =
          ZSTD_decompress(static_cast<void *>(buffer.data()), buffer.size(),
                          static_cast<const void *>(metadata.data()), metadata.size());
      if (ZSTD_isError(decompressed_size)) {  // NOLINT(*-implicit-bool-conversion)
        throw std::runtime_error(fmt::format("failed to decompress metadata using zstd: {}",
                                             ZSTD_getErrorName(decompressed_size)));
      }
      buffer.resize(decompressed_size);
      return buffer;
    }

    throw std::runtime_error(fmt::format("unknown compression method \"{}\"", compression));
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to read file metadata: {}", e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read file metadata: unknown error");
  }
}

[[nodiscard]] static glz::json_t parse_metadata_or_throw(const std::string &s) {
  try {
    return parse_json_string(s);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("not a valid JSON string: {}", e.what()));
  }
}

[[nodiscard]] static std::string import_metadata_from_parquet(
    const std::shared_ptr<arrow::io::ReadableFile> &fp,
    const std::vector<std::string> &ignored_keys) {
  try {
    auto json = parse_metadata_or_throw(read_metadata_or_throw(fp));
    json = strip_keys_from_json(json, ignored_keys);
    return to_json_string(json);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to read file metadata: {}", e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read file metadata: unknown error");
  }
}

[[nodiscard]] static std::uint8_t parse_and_validate_format_version(
    const std::filesystem::path &path, const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  std::uint8_t format_version = 1;
  try {
    format_version = parse_numeric<std::uint8_t>(
        read_attribute_or_throw(get_file_schema(fp), "NCHG:format-version"));
  } catch (const std::out_of_range &) {  // NOLINT
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to read format-version attribute: {}", e.what()));
  } catch (...) {
    throw std::runtime_error("failed to read format-version attribute: unknown error");
  }

  if (format_version < 2) {
    throw std::runtime_error(fmt::format(
        "unable to open file \"{}\": reading files with version={} is no longer supported", path,
        format_version));
  }

  return format_version;
}

[[nodiscard]] static std::shared_ptr<const hictk::Reference> import_chromosomes_from_parquet(
    const std::shared_ptr<arrow::io::ReadableFile> &fp) {
  try {
    const auto metadata = parse_json_string(import_metadata_from_parquet(fp, {}));
    const auto chrom_str = to_json_string(metadata.at("chromosomes"));
    return std::make_shared<hictk::Reference>(parse_json_string<hictk::Reference>(chrom_str));
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

  auto reader = parquet::ParquetFileReader::Open(std::move(fp), props);
  return std::make_shared<parquet::StreamReader>(std::move(reader));
}

[[nodiscard]] static auto validate_record_type(
    const std::filesystem::path &path, const std::shared_ptr<arrow::io::ReadableFile> &fp,
    ParquetStatsFileReader::RecordType expected_type,
    const std::shared_ptr<arrow::Schema> &schema = nullptr) -> ParquetStatsFileReader::RecordType {
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
    const auto col_names = (!!schema ? schema : get_file_schema(fp))->field_names();

    auto file_has_nchg_compute_records = [&] {
      return std::ranges::equal(expected_columns_nchg_compute, col_names);
    };

    auto file_has_nchg_filter_records = [&] {
      return std::ranges::equal(expected_columns_nchg_filter, col_names);
    };

    switch (expected_type) {
      using enum ParquetStatsFileReader::RecordType;
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

ParquetStatsFileReader::ParquetStatsFileReader(std::filesystem::path path,
                                               const std::shared_ptr<arrow::io::ReadableFile> &fp,
                                               RecordType record_type, std::size_t buffer_size)
    : ParquetStatsFileReader(std::move(path), fp, import_chromosomes_from_parquet(fp), record_type,
                             buffer_size) {}

ParquetStatsFileReader::ParquetStatsFileReader(std::filesystem::path path,
                                               std::shared_ptr<arrow::io::ReadableFile> fp,
                                               std::shared_ptr<const hictk::Reference> chromosomes,
                                               RecordType record_type, std::size_t buffer_size)
    : _path(std::move(path)),
      _schema(get_file_schema(fp)),
      _type(validate_record_type(_path, fp, record_type, _schema)),
      _chroms(std::move(chromosomes)),
      // It's important that we parse the format version before attempting to read the file metadata
      _format_version(parse_and_validate_format_version(_path, fp)),
      _metadata(import_metadata_from_parquet(fp, {})),
      _sr(init_parquet_stream_reader(std::move(fp), buffer_size)) {}

ParquetStatsFileReader::ParquetStatsFileReader(const std::filesystem::path &path,
                                               RecordType record_type, std::size_t buffer_size)
    : ParquetStatsFileReader(path, open_parquet_file(path), record_type, buffer_size) {}

const std::filesystem::path &ParquetStatsFileReader::path() const noexcept { return _path; }

std::shared_ptr<arrow::Schema> ParquetStatsFileReader::file_schema() const noexcept {
  return _schema;
}

auto ParquetStatsFileReader::record_type() const noexcept -> RecordType { return _type; }

std::string_view ParquetStatsFileReader::metadata() const noexcept { return _metadata; }

std::string ParquetStatsFileReader::read_metadata(const std::filesystem::path &path,
                                                  const std::vector<std::string> &ignored_keys) {
  const auto fp = open_parquet_file(path);
  std::ignore = parse_and_validate_format_version(path, fp);
  return import_metadata_from_parquet(fp, ignored_keys);
}

std::uint8_t ParquetStatsFileReader::format_version() const noexcept { return _format_version; }

std::shared_ptr<const hictk::Reference> ParquetStatsFileReader::chromosomes() const noexcept {
  return _chroms;
}

}  // namespace nchg
