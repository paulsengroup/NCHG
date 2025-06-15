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

#include <arrow/array/array_base.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <arrow/record_batch.h>
#include <arrow/util/key_value_metadata.h>
#include <parallel_hashmap/btree.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/reference.hpp>
#include <memory>
#include <optional>
#include <string_view>

namespace nchg {

class ParquetStatsFileWriter {
  std::filesystem::path _path;
  std::shared_ptr<parquet::WriterProperties> _props;
  std::shared_ptr<arrow::io::FileOutputStream> _fp;
  std::unique_ptr<parquet::arrow::FileWriter> _writer;
  std::shared_ptr<const arrow::KeyValueMetadata> _metadata;

  std::size_t _chunk_size{};
  std::size_t _chunk_capacity{};
  std::size_t _size{};

  hictk::Reference _chroms;

  arrow::StringDictionary32Builder _chrom1{};
  arrow::UInt32Builder _start1;
  arrow::UInt32Builder _end1;

  arrow::StringDictionary32Builder _chrom2{};
  arrow::UInt32Builder _start2;
  arrow::UInt32Builder _end2;

  arrow::DoubleBuilder _pvalue;
  std::optional<arrow::DoubleBuilder> _pvalue_corrected;
  arrow::UInt64Builder _observed;
  arrow::DoubleBuilder _expected;
  arrow::DoubleBuilder _log_ratio;
  arrow::DoubleBuilder _odds;
  arrow::DoubleBuilder _omega;

 public:
  static constexpr std::uint8_t format_version{2};
  ParquetStatsFileWriter() = delete;
  ParquetStatsFileWriter(hictk::Reference chroms, const std::filesystem::path &path, bool force,
                         std::string_view compression_method, std::uint8_t compression_lvl,
                         std::size_t threads, const std::string &metadata,
                         std::size_t batch_size = 1'000'000);
  ParquetStatsFileWriter(const ParquetStatsFileWriter &other) = delete;
  ParquetStatsFileWriter(ParquetStatsFileWriter &&other) noexcept = delete;

  ~ParquetStatsFileWriter() noexcept;

  ParquetStatsFileWriter &operator=(const ParquetStatsFileWriter &other) = delete;
  ParquetStatsFileWriter &operator=(ParquetStatsFileWriter &&other) noexcept = delete;

  template <typename Record>
  void append(const Record &r);

  // If finalize() is not called before a ParquetStatsFileWriter writer object goes out of scope
  // the underlying file will be removed by the destructor
  void finalize();

  // Same as the above function, but can be used to finalize an empty file
  template <typename Record>
  void finalize();

 private:
  template <typename Record>
  void initialize_writer();

  void write_chunk();
  void reset_builders();

  template <typename ArrayBuilder, typename T>
  void append(ArrayBuilder &builder, const T &data);

  [[nodiscard]] std::shared_ptr<arrow::RecordBatch> finalize_chunk();

  template <typename ArrayBuilder>
  [[nodiscard]] static std::shared_ptr<arrow::Array> finalize_chunk(ArrayBuilder &builder);

  static void setup_threadpool(std::size_t size);
  void init_chrom_writers();
};

}  // namespace nchg

#include "../../parquet_stats_file_writer_impl.hpp"
