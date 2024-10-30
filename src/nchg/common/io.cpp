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

// clang-format off
#include <nchg/nchg.hpp>


#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/record_batch.h>
#include <arrow/type.h>
#include <arrow/io/file.h>
#include <arrow/util/key_value_metadata.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/encoding.h>
#include <parquet/properties.h>
#include <parquet/stream_reader.h>
#include <parallel_hashmap/phmap.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include "nchg/tools/io.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <hictk/genomic_interval.hpp>
#include <hictk/numeric_utils.hpp>
#include <hictk/reference.hpp>
#include <iosfwd>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nchg {

static std::shared_ptr<arrow::Array> make_chrom_dict(const hictk::Reference &chroms) {
  arrow::StringBuilder builder{};
  for (const auto &chrom : chroms) {
    if (!chrom.is_all()) {
      const auto status = builder.Append(std::string{chrom.name()});
      if (!status.ok()) {
        throw std::runtime_error(status.ToString());
      }
    }
  }

  auto result = builder.Finish();
  if (!result.status().ok()) {
    throw std::runtime_error(result.status().ToString());
  }

  return result.MoveValueUnsafe();
}

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
  // clang-format off
  static constexpr std::array<std::string_view, 12> expected_columns_nchg_compute{
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
  };

  static constexpr std::array<std::string_view, 13> expected_columns_nchg_filter{
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
  };
  // clang-format on

  try {
    const auto schema = get_file_schema(fp);

    auto schema_has_column = [&](const auto &col) {
      return !!schema->GetFieldByName(std::string{col});
    };

    auto file_has_nchg_compute_records = [&] {
      return std::ranges::all_of(expected_columns_nchg_compute, schema_has_column);
    };

    auto file_has_nchg_filter_records = [&] {
      return std::ranges::all_of(expected_columns_nchg_compute, schema_has_column);
    };

    switch (expected_type) {
      using T = ParquetStatsFile::RecordType;
      case T::NCHGCompute: {
        if (file_has_nchg_compute_records()) {
          return T::NCHGCompute;
        }
        if (file_has_nchg_filter_records()) {
          throw std::runtime_error(
              "unexpected record type: expected NCHGCompute records, found NCHGFilter");
        }
        throw std::runtime_error("unexpected record type: file was likely not generated by NCHG");
      }
      case T::NCHGFilter: {
        if (file_has_nchg_filter_records()) {
          return T::NCHGFilter;
        }
        if (file_has_nchg_compute_records()) {
          throw std::runtime_error(
              "unexpected record type: expected NCHGFilter records, found NCHGCompute");
        }
        throw std::runtime_error("unexpected record type: file was likely not generated by NCHG");
      }
      case T::infer: {
        if (file_has_nchg_compute_records()) {
          return T::NCHGCompute;
        }
        if (file_has_nchg_filter_records()) {
          return T::NCHGFilter;
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

RecordBatchBuilder::RecordBatchBuilder(hictk::Reference chroms) : _chroms(std::move(chroms)) {
  const auto dict = make_chrom_dict(_chroms);
  auto status = _chrom1.InsertMemoValues(*dict);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }

  status = _chrom2.InsertMemoValues(*dict);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

std::size_t RecordBatchBuilder::size() const noexcept { return _i; }
std::size_t RecordBatchBuilder::capacity() const noexcept {
  return static_cast<std::size_t>(_chrom1.capacity());
}

void RecordBatchBuilder::reset() {
  _chrom1.Reset();
  _start1.Reset();
  _end1.Reset();

  _chrom2.Reset();
  _start2.Reset();
  _end2.Reset();

  _pvalue.Reset();
  _pvalue_corrected.Reset();
  _observed.Reset();
  _expected.Reset();
  _log_ratio.Reset();
  _odds.Reset();
  _omega.Reset();

  _i = 0;
}

std::shared_ptr<arrow::RecordBatch> RecordBatchBuilder::get() {
  std::vector<std::shared_ptr<arrow::Array>> columns{};
  columns.reserve(13);

  columns.emplace_back(finish(_chrom1));
  columns.emplace_back(finish(_start1));
  columns.emplace_back(finish(_end1));

  columns.emplace_back(finish(_chrom2));
  columns.emplace_back(finish(_start2));
  columns.emplace_back(finish(_end2));

  columns.emplace_back(finish(_pvalue));

  if (_pvalue_corrected.length() != 0) {
    columns.emplace_back(finish(_pvalue_corrected));
  }

  columns.emplace_back(finish(_observed));
  columns.emplace_back(finish(_expected));

  columns.emplace_back(finish(_log_ratio));
  columns.emplace_back(finish(_odds));
  columns.emplace_back(finish(_omega));

  if (columns.size() == 13) {
    return arrow::RecordBatch::Make(internal::get_schema_padj(_chroms),
                                    static_cast<std::int64_t>(size()), columns);
  }
  return arrow::RecordBatch::Make(internal::get_schema(_chroms), static_cast<std::int64_t>(size()),
                                  columns);
}

void RecordBatchBuilder::write(parquet::arrow::FileWriter &writer) {
  const auto batch = get();
  const auto status = writer.WriteRecordBatch(*batch);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
  reset();
}

std::string_view truncate_bed3_record(std::string_view record, char sep) {
  const auto pos1 = record.find(sep);
  if (pos1 == std::string_view::npos) {
    throw std::runtime_error("invalid bed record, expected 3 tokens, found 1");
  }
  const auto pos2 = record.find('\t', pos1 + 1);
  if (pos2 == std::string_view::npos) {
    throw std::runtime_error("invalid bed record, expected 3 tokens, found 2");
  }
  const auto pos3 = record.find('\t', pos2 + 1);

  return record.substr(0, pos3);
}

phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, const std::filesystem::path &path) {
  if (path.empty()) {
    return {};
  }

  SPDLOG_INFO("reading the user-provided bin mask from {}...", path);
  phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> mask{};
  std::string buffer{};

  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  std::size_t i = 1;
  try {
    fs.open(path);

    for (; std::getline(fs, buffer); ++i) {
      if (buffer.empty()) {
        continue;
      }

      if (buffer.back() == '\r') {
        buffer.resize(buffer.size() - 1);
      }

      try {
        const auto record = truncate_bed3_record(buffer);
        auto domain = hictk::GenomicInterval::parse_bed(chroms, record);

        const auto num_bins = (domain.chrom().size() + bin_size - 1) / bin_size;
        auto match = mask.try_emplace(domain.chrom(), std::vector<bool>(num_bins, false));

        const std::size_t j0 = domain.start() / bin_size;
        const std::size_t j1 = (domain.end() / bin_size) + 1;

        for (std::size_t j = j0; j < j1; ++j) {
          match.first->second[j] = true;
        }

      } catch (const std::exception &e) {
        throw std::runtime_error(
            fmt::format("found an invalid record at line {} of file {}: {}", i, path, e.what()));
      }
    }

  } catch (const std::exception &) {
    if (!fs.eof()) {
      throw;
    }
  }

  std::size_t num_bad_bins = 0;
  for (const auto &[_, v] : mask) {
    num_bad_bins += std::accumulate(v.begin(), v.end(), 0uz);
  }

  SPDLOG_INFO("masked {} bad bins based on {} domains read from {}...", num_bad_bins, i - 1, path);
  return mask;
}

namespace internal {
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
      arrow::field("observed_count", arrow::uint32(),  false),
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
      arrow::field("observed_count",   arrow::uint32(),  false),
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

}  // namespace internal

}  // namespace nchg
