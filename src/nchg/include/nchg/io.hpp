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
#include <arrow/util/key_value_metadata.h>
#include <arrow/util/thread_pool.h>
#include <fmt/format.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/encoding.h>
#include <parquet/stream_reader.h>
#include <spdlog/spdlog.h>

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

// https://stackoverflow.com/a/16000226
template <typename T, typename = int>
struct has_pval_corrected : std::false_type {};

template <typename T>
struct has_pval_corrected<T, decltype((void)T::pval_corrected, 0)> : std::true_type {};

template <typename T, typename = int>
struct has_log_ratio : std::false_type {};

template <typename T>
struct has_log_ratio<T, decltype((void)T::log_ratio, 0)> : std::true_type {};

inline std::shared_ptr<arrow::Array> make_chrom_dict(const hictk::Reference &chroms) {
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

template <typename Stats>
class ParquetStatsFile {
  std::shared_ptr<const hictk::Reference> _chroms{};
  std::shared_ptr<parquet::StreamReader> _sr{};

 public:
  class iterator;
  ParquetStatsFile() = default;
  explicit ParquetStatsFile(const std::filesystem::path &path) {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path));

    auto props = parquet::default_reader_properties();
    props.set_buffer_size(1'000'000);

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

    auto schema = batch_reader->schema();
    const auto metadata = schema->field(0)->metadata();
    const auto chrom_names = metadata->keys();
    std::vector<std::uint32_t> chrom_sizes(chrom_names.size(), 0);
    std::transform(metadata->values().begin(), metadata->values().end(), chrom_sizes.begin(),
                   [](const auto &size) {
                     return hictk::internal::parse_numeric_or_throw<std::uint32_t>(size);
                   });
    _chroms = std::make_shared<const hictk::Reference>(chrom_names.begin(), chrom_names.end(),
                                                       chrom_sizes.begin());

    _sr = std::make_shared<parquet::StreamReader>(parquet::ParquetFileReader::Open(fp, props));
  }

  [[nodiscard]] static bool is_nchg_compute_parquet(const std::filesystem::path &path) {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path));

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

    auto schema = batch_reader->schema();

    // clang-format off
    const std::vector<std::string> expected_columns{
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
    // clang-format on

    if (static_cast<std::size_t>(schema->num_fields()) != expected_columns.size()) {
      return false;
    }

    for (std::int64_t i = 0; i < schema->num_fields(); ++i) {
      if (!schema->GetFieldByName(expected_columns[static_cast<std::size_t>(i)])) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] static bool is_nchg_filter_parquet(const std::filesystem::path &path) {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path));

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

    auto schema = batch_reader->schema();

    // clang-format off
    const std::vector<std::string> expected_columns{
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

    if (static_cast<std::size_t>(schema->num_fields()) != expected_columns.size()) {
      return false;
    }

    for (std::int64_t i = 0; i < schema->num_fields(); ++i) {
      if (!schema->GetFieldByName(expected_columns[static_cast<std::size_t>(i)])) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] std::shared_ptr<const hictk::Reference> chromosomes() const noexcept {
    return _chroms;
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
      const auto chrom1 = !!_chroms ? _chroms->at(*_buffer) : hictk::Chromosome{0, *_buffer, 1};
      *_sr >> start1;
      *_sr >> end1;

      *_sr >> *_buffer;
      const auto chrom2 = !!_chroms ? _chroms->at(*_buffer) : hictk::Chromosome{0, *_buffer, 1};
      *_sr >> start2;
      *_sr >> end2;

      *_sr >> _value.pval;
      if constexpr (has_pval_corrected<Stats>()) {
        *_sr >> _value.pval_corrected;
      }
      *_sr >> observed_count;
      *_sr >> _value.expected;
      if constexpr (has_log_ratio<Stats>()) {
        *_sr >> _value.log_ratio;
      }
      *_sr >> _value.odds_ratio;
      *_sr >> _value.omega;
      *_sr >> parquet::EndRow;

      _value.pixel = hictk::Pixel{chrom1, start1, end1, chrom2, start2, end2, observed_count};

      return *this;
    }
  };
};

[[nodiscard]] inline parquet::Compression::type parse_parquet_compression(std::string_view method) {
  if (method == "zstd") {
    return parquet::Compression::ZSTD;
  }
  if (method == "lz4") {
    return parquet::Compression::LZ4;
  }
  throw std::runtime_error(
      fmt::format(FMT_STRING("unrecognized compression method \"{}\""), method));
}

[[nodiscard]] static std::shared_ptr<arrow::Schema> get_schema(const hictk::Reference &chroms) {
  auto chrom_dtype = arrow::dictionary(arrow::int32(), arrow::utf8());

  std::vector<std::string> keys{};
  std::vector<std::string> values{};
  for (const auto &chrom : chroms) {
    if (!chrom.is_all()) {
      keys.emplace_back(std::string{chrom.name()});
      values.emplace_back(fmt::to_string(chrom.size()));
    }
  }

  auto metadata = std::make_shared<arrow::KeyValueMetadata>(keys, values);

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

[[nodiscard]] static std::shared_ptr<arrow::Schema> get_schema_padj(
    const hictk::Reference &chroms) {
  auto chrom_dtype = arrow::dictionary(arrow::int32(), arrow::utf8());

  std::vector<std::string> keys{};
  std::vector<std::string> values{};
  for (const auto &chrom : chroms) {
    if (!chrom.is_all()) {
      keys.emplace_back(std::string{chrom.name()});
      values.emplace_back(fmt::to_string(chrom.size()));
    }
  }

  auto metadata = std::make_shared<arrow::KeyValueMetadata>(keys, values);

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

class RecordBatchBuilder {
  std::size_t _i{};

  hictk::Reference _chroms{};

  arrow::StringDictionary32Builder _chrom1{};
  arrow::UInt32Builder _start1{};
  arrow::UInt32Builder _end1{};

  arrow::StringDictionary32Builder _chrom2{};
  arrow::UInt32Builder _start2{};
  arrow::UInt32Builder _end2{};

  arrow::DoubleBuilder _pvalue{};
  arrow::DoubleBuilder _pvalue_corrected{};
  arrow::UInt32Builder _observed{};
  arrow::DoubleBuilder _expected{};
  arrow::DoubleBuilder _log_ratio{};
  arrow::DoubleBuilder _odds{};
  arrow::DoubleBuilder _omega{};

 public:
  RecordBatchBuilder(const hictk::Reference &chroms) : _chroms(chroms) {
    auto dict = make_chrom_dict(_chroms);
    auto status = _chrom1.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    status = _chrom2.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  }

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
    append(_log_ratio, s.log_ratio);
    append(_odds, s.odds_ratio);
    append(_omega, s.omega);

    if constexpr (has_pval_corrected<Stats>::value) {
      append(_pvalue_corrected, s.pval_corrected);
    }

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
    _pvalue_corrected.Reset();
    _observed.Reset();
    _expected.Reset();
    _log_ratio.Reset();
    _odds.Reset();
    _omega.Reset();

    _i = 0;
  }

  [[nodiscard]] std::shared_ptr<arrow::RecordBatch> get() {
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
      return arrow::RecordBatch::Make(get_schema_padj(_chroms), static_cast<std::int64_t>(size()),
                                      columns);
    }
    return arrow::RecordBatch::Make(get_schema(_chroms), static_cast<std::int64_t>(size()),
                                    columns);
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

template <typename Record>
[[nodiscard]] static std::unique_ptr<parquet::arrow::FileWriter> init_parquet_file_writer(
    const hictk::Reference &chroms, const std::filesystem::path &path, bool force,
    std::string_view compression_method, std::uint8_t compression_lvl, std::size_t threads) {
  if (path.empty()) {
    return {};
  }

  const auto schema =
      has_pval_corrected<Record>::value ? *get_schema_padj(chroms) : *get_schema(chroms);

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

  auto arrow_properties = parquet::ArrowWriterProperties::Builder()
                              .set_use_threads(threads > 1)
                              ->store_schema()
                              ->build();

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
