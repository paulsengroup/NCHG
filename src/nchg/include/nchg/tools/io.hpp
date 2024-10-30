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
#include <arrow/builder.h>
#include <arrow/record_batch.h>
#include <arrow/type_fwd.h>
#include <parquet/stream_reader.h>
#include <parallel_hashmap/phmap.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/reference.hpp>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>

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

class ParquetStatsFile {
 public:
  enum class RecordType : std::uint_fast8_t { infer, NCHGCompute, NCHGFilter };
  template <typename Stats>
  class iterator;

 private:
  RecordType _type{RecordType::NCHGCompute};
  std::shared_ptr<const hictk::Reference> _chroms{};
  std::shared_ptr<parquet::StreamReader> _sr{};

  ParquetStatsFile(const std::filesystem::path &path,
                   const std::shared_ptr<arrow::io::ReadableFile> &fp, RecordType record_type,
                   std::size_t buffer_size);
  ParquetStatsFile(const std::filesystem::path &path, std::shared_ptr<arrow::io::ReadableFile> fp,
                   std::shared_ptr<const hictk::Reference> chromosomes, RecordType record_type,
                   std::size_t buffer_size);

 public:
  ParquetStatsFile() = default;
  ParquetStatsFile(const std::filesystem::path &path, RecordType record_type,
                   std::size_t buffer_size = 1'000'000);

  [[nodiscard]] auto record_type() const noexcept -> RecordType;

  [[nodiscard]] std::shared_ptr<const hictk::Reference> chromosomes() const noexcept;
  template <typename Stats>
  [[nodiscard]] auto begin() -> iterator<Stats>;
  template <typename Stats>
  [[nodiscard]] auto end() -> iterator<Stats>;

  template <typename Stats>
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

    iterator(std::shared_ptr<const hictk::Reference> chroms,
             std::shared_ptr<parquet::StreamReader> sr, bool init_value = true);
    [[nodiscard]] static auto at_end(std::shared_ptr<const hictk::Reference> chroms,
                                     std::shared_ptr<parquet::StreamReader> sr) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const noexcept -> const_reference;
    [[nodiscard]] auto operator->() const noexcept -> const_pointer;

    auto operator++() -> iterator &;

   private:
    void read_pixel();
  };
};

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
  RecordBatchBuilder(hictk::Reference chroms);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t capacity() const noexcept;

  template <typename Stats>
  void append(const Stats &s);
  void reset();

  [[nodiscard]] std::shared_ptr<arrow::RecordBatch> get();

  void write(parquet::arrow::FileWriter &writer);

 private:
  template <typename ArrayBuilder, typename T>
  void append(ArrayBuilder &builder, const T &data);

  template <typename ArrayBuilder>
  [[nodiscard]] std::shared_ptr<arrow::Array> finish(ArrayBuilder &builder);
};

template <typename Record>
[[nodiscard]] std::unique_ptr<parquet::arrow::FileWriter> init_parquet_file_writer(
    const hictk::Reference &chroms, const std::filesystem::path &path, bool force,
    std::string_view compression_method, std::uint8_t compression_lvl, std::size_t threads);

[[nodiscard]] phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, const std::filesystem::path &path);

[[nodiscard]] std::string_view truncate_bed3_record(std::string_view record, char sep = '\t');

namespace internal {

[[nodiscard]] std::shared_ptr<arrow::Schema> get_schema(const hictk::Reference &chroms);
[[nodiscard]] std::shared_ptr<arrow::Schema> get_schema_padj(const hictk::Reference &chroms);

[[nodiscard]] parquet::Compression::type parse_parquet_compression(std::string_view method);

}  // namespace internal

}  // namespace nchg

#include "../../../common/io_impl.hpp"
