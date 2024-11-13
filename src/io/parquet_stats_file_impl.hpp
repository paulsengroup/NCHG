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
#include <arrow/array.h>
#include <arrow/io/file.h>
#include <arrow/util/thread_pool.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "nchg/parquet_helpers.hpp"

namespace nchg {

template <typename Stats>
inline auto ParquetStatsFile::begin() -> iterator<Stats> {
  return {_chroms, _sr, true};
}
template <typename Stats>
inline auto ParquetStatsFile::end() -> iterator<Stats> {
  return iterator<Stats>::at_end(_chroms, _sr);
}

template <typename Stats>
inline ParquetStatsFile::iterator<Stats>::iterator(std::shared_ptr<const hictk::Reference> chroms,
                                                   std::shared_ptr<parquet::StreamReader> sr,
                                                   bool init_value)
    : _chroms(std::move(chroms)), _sr(std::move(sr)), _buffer(std::make_shared<std::string>()) {
  if (init_value && _sr->current_row() != _sr->num_rows()) {
    if (_sr->eof()) {
      *this = at_end(_chroms, _sr);
    } else {
      read_pixel();
    }
  }
}

template <typename Stats>
inline auto ParquetStatsFile::iterator<Stats>::at_end(
    std::shared_ptr<const hictk::Reference> chroms, std::shared_ptr<parquet::StreamReader> sr)
    -> iterator<Stats> {
  iterator it{std::move(chroms), std::move(sr), false};
  it._offset = it._sr->num_rows();

  return it;
}

template <typename Stats>
inline bool ParquetStatsFile::iterator<Stats>::operator==(const iterator &other) const noexcept {
  return _sr == other._sr && _offset == other._offset;
}

template <typename Stats>
inline bool ParquetStatsFile::iterator<Stats>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename Stats>
inline auto ParquetStatsFile::iterator<Stats>::operator*() const noexcept -> const_reference {
  return _value;
}

template <typename Stats>
inline auto ParquetStatsFile::iterator<Stats>::operator->() const noexcept -> const_pointer {
  return &_value;
}

template <typename Stats>
inline auto ParquetStatsFile::iterator<Stats>::operator++() -> iterator & {
  if (_sr->eof()) [[unlikely]] {
    *this = at_end(_chroms, _sr);
    return *this;
  }

  read_pixel();
  return *this;
}

template <typename Stats>
inline void ParquetStatsFile::iterator<Stats>::read_pixel() {
  assert(!_sr->eof());
  std::uint32_t start1{};
  std::uint32_t end1{};
  std::uint32_t start2{};
  std::uint32_t end2{};
  std::uint64_t observed_count{};

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
  *_sr >> _value.log_ratio;

  *_sr >> _value.odds_ratio;
  *_sr >> _value.omega;
  *_sr >> parquet::EndRow;

  _value.pixel = hictk::Pixel{chrom1, start1, end1, chrom2, start2, end2, observed_count};
}

}  // namespace nchg
