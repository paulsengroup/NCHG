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
#include <arrow/util/thread_pool.h>
#include <parquet/properties.h>
#include <parquet/arrow/writer.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdint>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>

#include "nchg/parquet_helpers.hpp"

namespace nchg {

template <typename Record>
inline void ParquetStatsFileWriter::initialize_writer() {
  SPDLOG_DEBUG("initializing file \"{}\"...", _path.string());
  try {
    if (_writer) {
      throw std::logic_error("ParquetStatsFileWriter has already been initialized");
    }

    const auto schema =
        has_pval_corrected<Record>::value ? *get_schema_padj(_chroms) : *get_schema(_chroms);

    if constexpr (has_pval_corrected<Record>::value) {
      _pvalue_corrected = {};
    }

    auto arrow_properties = parquet::ArrowWriterProperties::Builder()
                                .set_use_threads(arrow::GetCpuThreadPoolCapacity() > 1)
                                ->store_schema()
                                ->build();

    auto result = parquet::arrow::FileWriter::Open(schema, arrow::default_memory_pool(), _fp,
                                                   _props, arrow_properties);

    if (!result.status().ok()) {
      throw std::runtime_error(result.status().ToString());
    }

    _writer = result.MoveValueUnsafe();
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to initialize the file writer for file \"{}\": {}", _path, e.what()));
  }
}

template <typename Record>
inline void ParquetStatsFileWriter::append(const Record &r) {
  if (!_writer) [[unlikely]] {
    initialize_writer<Record>();
  }

  if (_chunk_size == _chunk_capacity) {
    write_chunk();
  }

  append(_chrom1, std::string{r.pixel.coords.bin1.chrom().name()});
  append(_start1, r.pixel.coords.bin1.start());
  append(_end1, r.pixel.coords.bin1.end());

  append(_chrom2, std::string{r.pixel.coords.bin2.chrom().name()});
  append(_start2, r.pixel.coords.bin2.start());
  append(_end2, r.pixel.coords.bin2.end());

  append(_pvalue, r.pval);
  append(_observed, r.pixel.count);
  append(_expected, r.expected);
  append(_log_ratio, r.log_ratio);
  append(_odds, r.odds_ratio);
  append(_omega, r.omega);

  if constexpr (has_pval_corrected<Record>::value) {
    assert(_pvalue_corrected);
    append(*_pvalue_corrected, r.pval_corrected);
  } else {
    assert(!_pvalue_corrected);
  }
  ++_chunk_size;
}

template <typename Record>
inline void ParquetStatsFileWriter::finalize() {
  if (!_writer) {
    initialize_writer<Record>();
  }
  finalize();
}

template <typename ArrayBuilder, typename T>
inline void ParquetStatsFileWriter::append(ArrayBuilder &builder, const T &data) {
  const auto status = builder.Append(data);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

template <typename ArrayBuilder>
inline std::shared_ptr<arrow::Array> ParquetStatsFileWriter::finalize_chunk(ArrayBuilder &builder) {
  auto result = builder.Finish();
  if (!result.status().ok()) {
    throw std::runtime_error(result.status().ToString());
  }

  return result.MoveValueUnsafe();
}

}  // namespace nchg
