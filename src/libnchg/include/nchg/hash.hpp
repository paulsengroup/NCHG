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

#include <cstddef>
#include <functional>

namespace nchg::internal {
// Adapted from:
// https://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine

template <typename T>
[[nodiscard]] inline std::size_t hash_combine(std::size_t seed, const T &v) {
  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
  return seed;
}
template <typename T, typename... Args>
[[nodiscard]] inline std::size_t hash_combine(std::size_t seed, const T &v, const Args &...args) {
  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
  return hash_combine(seed, args...);
}
}  // namespace nchg::internal
