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

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <version>

namespace nchg {
// to avoid useless casts (see
// https://github.com/nlohmann/json/issues/2893#issuecomment-889152324)
template <typename T, typename U>
[[maybe_unused]] [[nodiscard]] constexpr T conditional_static_cast(U value) {
  if constexpr (std::is_same_v<T, U>) {
    return value;
  } else {
    return static_cast<T>(value);
  }
}

struct identity {
  template <typename T>
  [[nodiscard]] constexpr T &&operator()(T &&a) const noexcept {
    return std::forward<T>(a);
  }
  using is_transparent = void;
};

[[nodiscard]] constexpr bool ndebug_defined() noexcept {
#ifdef NDEBUG
  return true;
#else
  return false;
#endif
}

[[nodiscard]] constexpr bool ndebug_not_defined() noexcept { return !ndebug_defined(); }

#ifdef __cpp_lib_unreachable
#define NCHG_UNREACHABLE_CODE std::unreachable();
#elifdef __GNUC__
#define NCHG_UNREACHABLE_CODE __builtin_unreachable();
#elifdef _MSC_VER
#define NCHG_UNREACHABLE_CODE __assume(0);
#else
#define NCHG_UNREACHABLE_CODE std::abort();
#endif

[[noreturn]] inline void unreachable_code() {
  if constexpr (ndebug_not_defined()) {
    throw std::logic_error("Unreachable code");
  }
  NCHG_UNREACHABLE_CODE
}

}  // namespace nchg
