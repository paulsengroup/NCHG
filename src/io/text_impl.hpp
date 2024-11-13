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

#include <fmt/format.h>

#include <cstddef>
#include <stdexcept>
#include <string_view>

namespace nchg {

template <std::size_t NTOKS>
constexpr std::string_view truncate_record(std::string_view record, char sep) {
  static_assert(NTOKS != 0);

  std::size_t offset{};
  for (std::size_t i = 0; i < NTOKS; ++i) {
    const auto pos = record.find(sep, offset + 1);
    if (pos == std::string_view::npos && i != NTOKS - 1) [[unlikely]] {
      throw std::runtime_error(
          fmt::format("invalid record, expected {} tokens, found {}", NTOKS, i));
    }
    offset = pos;
  }

  return record.substr(0, offset);
}

}  // namespace nchg
