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

#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/io/concurrency.h>
#include <arrow/io/file.h>
#include <arrow/record_batch.h>
#include <arrow/util/thread_pool.h>
#include <fmt/format.h>
#include <moodycamel/blockingconcurrentqueue.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_writer.h>
#include <parquet/stream_reader.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "nchg/k_merger.hpp"

namespace nchg {
template <typename Stats>
class ParquetStatsFile {
  std::shared_ptr<const hictk::Reference> _chroms{};
  std::shared_ptr<parquet::StreamReader> _sr{};

 public:
  class iterator;

  ParquetStatsFile(hictk::Reference chroms, const std::filesystem::path &path)
      : _chroms(std::make_shared<const hictk::Reference>(std::move(chroms))) {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path));

    auto props = parquet::default_reader_properties();
    props.set_buffer_size(1'000'000);

    _sr = std::make_shared<parquet::StreamReader>(parquet::ParquetFileReader::Open(fp, props));
  }

  [[nodiscard]] auto begin() -> iterator { return {_chroms, _sr, true}; }
  [[nodiscard]] auto end() -> iterator { return iterator::at_end(_chroms, _sr); }

  class iterator {
    std::shared_ptr<const hictk::Reference> _chroms{};
    std::shared_ptr<parquet::StreamReader> _sr{};
    std::shared_ptr<std::string> _buffer{};
    Stats _value{};
    std::int64_t _offset{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Stats;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator(std::shared_ptr<const hictk::Reference> chroms,
             std::shared_ptr<parquet::StreamReader> sr, bool init_value = true)
        : _chroms(std::move(chroms)), _sr(std::move(sr)), _buffer(std::make_shared<std::string>()) {
      if (init_value && _sr->current_row() != _sr->num_rows()) {
        *this = ++*this;
      }
    }

    [[nodiscard]] static iterator at_end(std::shared_ptr<const hictk::Reference> chroms,
                                         std::shared_ptr<parquet::StreamReader> sr) {
      iterator it{std::move(chroms), std::move(sr), false};
      it._offset = it._sr->num_rows();

      return it;
    }

    [[nodiscard]] bool operator==(const iterator &other) const noexcept {
      return _sr == other._sr && _offset == other._offset;
    }
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept {
      return !(*this == other);
    }

    auto operator*() const noexcept -> const_reference { return _value; }
    auto operator->() const noexcept -> const_pointer { return &_value; }

    [[nodiscard]] auto operator++() -> iterator & {
      if (_sr->eof()) {
        *this = at_end(_chroms, _sr);
        return *this;
      }

      std::uint32_t start1{};
      std::uint32_t end1{};
      std::uint32_t start2{};
      std::uint32_t end2{};
      std::uint32_t observed_count{};

      *_sr >> *_buffer;
      const auto &chrom1 = _chroms->at(*_buffer);
      *_sr >> start1;
      *_sr >> end1;

      *_sr >> *_buffer;
      const auto &chrom2 = _chroms->at(*_buffer);
      *_sr >> start2;
      *_sr >> end2;

      *_sr >> _value.pval;
      *_sr >> observed_count;
      *_sr >> _value.expected;
      *_sr >> _value.odds_ratio;
      *_sr >> _value.omega;
      *_sr >> parquet::EndRow;

      _value.pixel = hictk::Pixel{chrom1, start1, end1, chrom2, start2, end2, observed_count};

      return *this;
    }
  };
};

[[nodiscard]] static parquet::Compression::type parse_parquet_compression(std::string_view method) {
  if (method == "zstd") {
    return parquet::Compression::ZSTD;
  }
  if (method == "lz4") {
    return parquet::Compression::LZ4;
  }
  throw std::runtime_error(
      fmt::format(FMT_STRING("unrecognized compression method \"{}\""), method));
}

[[nodiscard]] static std::shared_ptr<arrow::Schema> get_schema() {
  return arrow::schema({
      // clang-format off
      arrow::field("chrom1",         arrow::utf8()),
      arrow::field("start1",         arrow::uint32()),
      arrow::field("end1",           arrow::uint32()),
      arrow::field("chrom2",         arrow::utf8()),
      arrow::field("start2",         arrow::uint32()),
      arrow::field("end2",           arrow::uint32()),
      arrow::field("pvalue",         arrow::float64()),
      arrow::field("observed_count", arrow::uint32()),
      arrow::field("expected_count", arrow::float64()),
      arrow::field("odds_ratio",     arrow::float64()),
      arrow::field("omega",          arrow::float64())
      // clang-format on
  });
}

class RecordBatchBuilder {
  std::size_t _i{};

  arrow::StringBuilder _chrom1{};
  arrow::UInt32Builder _start1{};
  arrow::UInt32Builder _end1{};

  arrow::StringBuilder _chrom2{};
  arrow::UInt32Builder _start2{};
  arrow::UInt32Builder _end2{};

  arrow::DoubleBuilder _pvalue{};
  arrow::UInt32Builder _observed{};
  arrow::DoubleBuilder _expected{};
  arrow::DoubleBuilder _odds{};
  arrow::DoubleBuilder _omega{};

 public:
  RecordBatchBuilder() = default;

  [[nodiscard]] std::size_t size() const noexcept { return _i; }
  [[nodiscard]] std::size_t capacity() const noexcept {
    return static_cast<std::size_t>(_chrom1.capacity());
  }

  template <typename Stats>
  void append(const Stats &s) {
    append(_chrom1, std::string{s.pixel.coords.bin1.chrom().name()});
    append(_start1, s.pixel.coords.bin1.start());
    append(_end1, s.pixel.coords.bin1.end());

    append(_chrom2, std::string{s.pixel.coords.bin2.chrom().name()});
    append(_start2, s.pixel.coords.bin2.start());
    append(_end2, s.pixel.coords.bin2.end());

    append(_pvalue, s.pval);
    append(_observed, s.pixel.count);
    append(_expected, s.expected);
    append(_odds, s.odds_ratio);
    append(_omega, s.omega);

    ++_i;
  }

  void reset() {
    _chrom1.Reset();
    _start1.Reset();
    _end1.Reset();

    _chrom2.Reset();
    _start2.Reset();
    _end2.Reset();

    _pvalue.Reset();
    _observed.Reset();
    _expected.Reset();
    _odds.Reset();
    _omega.Reset();

    _i = 0;
  }

  [[nodiscard]] std::shared_ptr<arrow::RecordBatch> get() {
    return arrow::RecordBatch::Make(
        get_schema(), static_cast<std::int64_t>(size()),
        {
            // clang-format off
        finish(_chrom1),
        finish(_start1),
        finish(_end1),
        finish(_chrom2),
        finish(_start2),
        finish(_end2),
        finish(_pvalue),
        finish(_observed),
        finish(_expected),
        finish(_odds),
        finish(_omega)
            // clang-format on
        });
  }

  void write(parquet::arrow::FileWriter &writer) {
    const auto batch = get();
    const auto status = writer.WriteRecordBatch(*batch);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
    reset();
  }

 private:
  template <typename ArrayBuilder, typename T>
  void append(ArrayBuilder &builder, const T &data) {
    const auto status = builder.Append(data);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  }

  template <typename ArrayBuilder>
  [[nodiscard]] std::shared_ptr<arrow::Array> finish(ArrayBuilder &builder) {
    auto result = builder.Finish();
    if (!result.status().ok()) {
      throw std::runtime_error(result.status().ToString());
    }

    return result.MoveValueUnsafe();
  }
};

[[nodiscard]] static std::unique_ptr<parquet::arrow::FileWriter> init_parquet_file_writer(
    const std::filesystem::path &path, bool force, std::string_view compression_method,
    std::uint8_t compression_lvl, std::size_t threads) {
  if (path.empty()) {
    return {};
  }

  const auto schema = *get_schema();

  auto builder = parquet::WriterProperties::Builder()
                     .created_by("NCHG v0.0.1")
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

  auto arrow_properties =
      parquet::ArrowWriterProperties::Builder().set_use_threads(threads > 1)->build();

  if (force) {
    std::filesystem::remove(path);
  }

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
