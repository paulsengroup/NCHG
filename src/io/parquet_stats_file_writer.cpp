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

#include "nchg/parquet_stats_file_writer.hpp"

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/array/array_base.h>
#include <arrow/builder.h>
#include <arrow/record_batch.h>
#include <arrow/util/thread_pool.h>
#include <parquet/properties.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "nchg/parquet_stats_file_reader.hpp"

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

[[nodiscard]] static std::shared_ptr<arrow::io::FileOutputStream> create_parquet_file(
    const std::filesystem::path &path, bool force) {
  if (path.empty()) {
    return {};
  }

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
  return f;
}

ParquetStatsFileWriter::ParquetStatsFileWriter(hictk::Reference chroms,
                                               const std::filesystem::path &path, bool force,
                                               std::string_view compression_method,
                                               std::uint8_t compression_lvl, std::size_t threads,
                                               std::size_t batch_size)
    : _path(path),
      _props(parquet::WriterProperties::Builder()
                 .created_by("NCHG v0.0.2")  // TODO fixme
                 ->version(parquet::ParquetVersion::PARQUET_2_6)
                 ->data_page_version(parquet::ParquetDataPageVersion::V2)
                 ->compression(parse_parquet_compression(compression_method))
                 ->compression_level(compression_lvl)
                 ->disable_statistics()
                 ->build()),
      _fp(create_parquet_file(path, force)),
      _chunk_capacity(batch_size),
      _chroms(std::move(chroms)) {
  if (_chunk_capacity == 0) {
    throw std::logic_error("batch_size cannot be 0");
  }

  setup_threadpool(threads);
  init_chrom_writers();
}

ParquetStatsFileWriter::~ParquetStatsFileWriter() noexcept {
  if (!_fp) {
    return;
  }

  try {
    SPDLOG_DEBUG("removing file \"{}\" because it was never finalized...", _path.string());
    _writer.reset();
    _fp.reset();
    std::filesystem::remove(_path);  // NOLINT
  } catch (const std::exception &e) {
    SPDLOG_ERROR("failed to remove file \"{}\": {}", _path.string(), e.what());
  } catch (...) {
    SPDLOG_ERROR("failed to remove file \"{}\": unknown error", _path.string());
  }
}

void ParquetStatsFileWriter::finalize() {
  if (_writer) {
    if (_chunk_size != 0) {
      write_chunk();
    }
    _writer.reset();
    _fp.reset();
    return;
  }

  if (_size == 0) {
    throw std::logic_error(
        "finalize() was called on a ParquetStatsFileWriter instance that was never initialized!");
  }

  throw std::logic_error(
      "finalize() was called on a ParquetStatsFileWriter instance that has already been "
      "finalized!");
}

void ParquetStatsFileWriter::setup_threadpool(std::size_t size) {
  try {
    if (size > 1 && static_cast<std::int32_t>(size - 1) > arrow::GetCpuThreadPoolCapacity()) {
      auto status = arrow::SetCpuThreadPoolCapacity(static_cast<std::int32_t>(size - 1));
      if (!status.ok()) {
        throw std::runtime_error(status.ToString());
      }
      status = arrow::io::SetIOThreadPoolCapacity(static_cast<std::int32_t>(size - 1));
      if (!status.ok()) {
        throw std::runtime_error(status.ToString());
      }
    }
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(fmt::format("failed to setup Arrow's thread pool: {}", e.what()));
  }
}

void ParquetStatsFileWriter::init_chrom_writers() {
  try {
    const auto dict = make_chrom_dict(_chroms);
    auto status = _chrom1.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    status = _chrom2.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(
        fmt::format("failed to setup one of the chromosome builders: {}", e.what()));
  }
}

std::shared_ptr<arrow::RecordBatch> ParquetStatsFileWriter::finalize_chunk() {
  std::vector<std::shared_ptr<arrow::Array>> columns(13);
  columns.clear();

  columns.emplace_back(finalize_chunk(_chrom1));
  columns.emplace_back(finalize_chunk(_start1));
  columns.emplace_back(finalize_chunk(_end1));

  columns.emplace_back(finalize_chunk(_chrom2));
  columns.emplace_back(finalize_chunk(_start2));
  columns.emplace_back(finalize_chunk(_end2));

  columns.emplace_back(finalize_chunk(_pvalue));

  if (_pvalue_corrected) {
    columns.emplace_back(finalize_chunk(*_pvalue_corrected));
  }

  columns.emplace_back(finalize_chunk(_observed));
  columns.emplace_back(finalize_chunk(_expected));

  columns.emplace_back(finalize_chunk(_log_ratio));
  columns.emplace_back(finalize_chunk(_odds));
  columns.emplace_back(finalize_chunk(_omega));

  if (_pvalue_corrected) {
    return arrow::RecordBatch::Make(get_schema_padj(_chroms),
                                    static_cast<std::int64_t>(_chunk_size), columns);
  }
  return arrow::RecordBatch::Make(get_schema(_chroms), static_cast<std::int64_t>(_chunk_size),
                                  columns);
}

void ParquetStatsFileWriter::write_chunk() {
  try {
    SPDLOG_DEBUG("writing chunk #{} to file \"{}\"...", (_size / _chunk_size) + 1, _path.string());
    assert(_writer);
    const auto batch = finalize_chunk();
    const auto status = _writer->WriteRecordBatch(*batch);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
    _size += _chunk_size;
    reset_builders();
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to write chunk to file \"{}\": {}", _path.string(), e.what()));
  }
}

void ParquetStatsFileWriter::reset_builders() {
  _chrom1.Reset();
  _start1.Reset();
  _end1.Reset();

  _chrom2.Reset();
  _start2.Reset();
  _end2.Reset();

  _pvalue.Reset();

  if (_pvalue_corrected) {
    _pvalue_corrected->Reset();
  }
  _observed.Reset();
  _expected.Reset();
  _log_ratio.Reset();
  _odds.Reset();
  _omega.Reset();

  _chunk_size = 0;
}

}  // namespace nchg
