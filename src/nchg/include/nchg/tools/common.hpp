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

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <chrono>
#include <string>

namespace nchg {
namespace internal {
[[nodiscard]] constexpr std::string strip_leading_zero(std::string s) {
  if (s.front() != '0') {
    return s;
  }
  return s.substr(1);
}
}  // namespace internal

template <typename Duration>
[[nodiscard]] constexpr std::string format_duration(const Duration& duration) {
  if (duration < std::chrono::microseconds(1)) {
    return fmt::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(duration));
  }

  if (duration < std::chrono::milliseconds(1)) {
    return fmt::to_string(std::chrono::duration_cast<std::chrono::microseconds>(duration));
  }

  if (duration < std::chrono::seconds(1)) {
    return fmt::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(duration));
  }

  if (duration < std::chrono::minutes(1)) {
    return internal::strip_leading_zero(
        fmt::format("{:%S}s", std::chrono::duration_cast<std::chrono::milliseconds>(duration)));
  }

  if (duration < std::chrono::hours(1)) {
    return internal::strip_leading_zero(
        fmt::format("{:%Mm:%S}s", std::chrono::duration_cast<std::chrono::milliseconds>(duration)));
  }

  if (duration < std::chrono::days(1)) {
    return internal::strip_leading_zero(fmt::format(
        "{:%Hh:%Mm:%S}s", std::chrono::duration_cast<std::chrono::milliseconds>(duration)));
  }

  return fmt::format("{:%jd:%Hh:%Mm:%S}s",
                     std::chrono::duration_cast<std::chrono::milliseconds>(duration));
}

}  // namespace nchg
