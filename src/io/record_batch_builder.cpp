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

#include "nchg/record_batch_builder.hpp"

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/array/array_base.h>
#include <arrow/builder.h>
#include <arrow/record_batch.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <cstddef>
#include <cstdint>
#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "nchg/parquet_stats_file.hpp"

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
    return arrow::RecordBatch::Make(get_schema_padj(_chroms), static_cast<std::int64_t>(size()),
                                    columns);
  }
  return arrow::RecordBatch::Make(get_schema(_chroms), static_cast<std::int64_t>(size()), columns);
}

void RecordBatchBuilder::write(parquet::arrow::FileWriter &writer) {
  const auto batch = get();
  const auto status = writer.WriteRecordBatch(*batch);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
  reset();
}

}  // namespace nchg
