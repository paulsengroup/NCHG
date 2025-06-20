// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
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

#include <string>
#include <string_view>

namespace glz {
struct json_t;
}

namespace nchg {

template <typename T>
[[nodiscard]] std::string to_json_string(const T &value);

template <typename T>
void to_json_string(const T &value, std::string &buff);

template <typename T>
[[nodiscard]] T parse_json_string(std::string_view buff);

[[nodiscard]] glz::json_t parse_json_string(std::string_view buff);

}  // namespace nchg

#include "../../metadata_impl.hpp"
