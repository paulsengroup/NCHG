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
#include <arrow/builder.h>
#include <arrow/record_batch.h>
#include <arrow/type_fwd.h>
#include <parquet/arrow/writer.h>
#include <parallel_hashmap/phmap.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/reference.hpp>
#include <memory>
#include <string_view>
#include <type_traits>

namespace nchg {

template <typename Record>
[[nodiscard]] std::unique_ptr<parquet::arrow::FileWriter> init_parquet_file_writer(
    const hictk::Reference &chroms, const std::filesystem::path &path, bool force,
    std::string_view compression_method, std::uint8_t compression_lvl, std::size_t threads);

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
  arrow::UInt64Builder _observed{};
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

}  // namespace nchg

#include "../../record_batch_builder_impl.hpp"
