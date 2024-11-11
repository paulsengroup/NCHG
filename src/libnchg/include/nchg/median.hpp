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

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace nchg {

template <typename T>
[[nodiscard]] constexpr T median(std::vector<T> v) {
  if (v.empty()) {
    throw std::runtime_error("median was called on an empty vector");
  }

  const auto size = static_cast<std::ptrdiff_t>(v.size());
  auto first = v.begin();
  auto mid = first + (size / 2);
  auto last = v.end();

  std::nth_element(first, mid, last);

  if (size % 2 != 0) {
    return *mid;
  }

  const auto n1 = *mid;
  std::nth_element(first, --mid, last);
  const auto n2 = *mid;

  return (n1 + n2) / 2;
}

}  // namespace nchg
