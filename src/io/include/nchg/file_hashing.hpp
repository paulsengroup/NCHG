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
NCHG_DISABLE_WARNING_MACRO_REDEFINED
NCHG_DISABLE_WARNING_OLD_STYLE_CAST
#include <xxh3.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <compare>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

namespace nchg {

class XXH3Digest {
  std::string _digest;

 public:
  static constexpr std::size_t WIDTH = 128;

  struct glaze {
    friend class XXH3Digest;
    using T = XXH3Digest;
    static constexpr auto value = &T::_digest;
  };

  XXH3Digest();
  explicit XXH3Digest(std::string digest);

  [[nodiscard]] std::string_view operator()() const noexcept;

  std::strong_ordering operator<=>(const XXH3Digest& other) const noexcept = default;
};

struct XXH3StateDeleter {
  void operator()(XXH3_state_t* state) const noexcept;
};

using XXH3StatePtr = std::unique_ptr<XXH3_state_t, XXH3StateDeleter>;

[[nodiscard]] XXH3StatePtr init_xxh3_state();

void update_state(XXH3StatePtr& state, std::string_view buffer);

[[nodiscard]] std::string to_hex(const XXH3StatePtr& state);

std::string hash_file(const std::filesystem::path& path, std::streamsize sample_size);

}  // namespace nchg
