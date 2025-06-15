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

#include <algorithm>
#include <glaze/glaze.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/reference.hpp>
#include <iterator>
#include <stdexcept>
#include <string>

namespace glz {

namespace internal::nchg {
struct Chromosome {
  std::uint32_t id;
  std::string name;
  std::uint32_t size;
};

}  // namespace internal::nchg

template <>
struct meta<internal::nchg::Chromosome> {
  using T = internal::nchg::Chromosome;
  static constexpr auto value = object(
      // clang-format off
      "id", &T::id,
      "name", &T::name,
      "size", &T::size
      // clang-format on
  );
};

template <>
struct meta<hictk::Chromosome> {
  using T = hictk::Chromosome;
  static constexpr auto value = object(
      // clang-format off
      "id", [](const T& chrom) { return chrom.id(); },
      "name", [](const T& chrom) { return chrom.name(); },
      "size", [](const T& chrom) { return chrom.size(); }
      // clang-format on
  );
};

template <>
struct from<JSON, hictk::Reference> {
  template <auto Opts>
  static void op(hictk::Reference& chroms, auto&&... args) {
    std::vector<internal::nchg::Chromosome> buff;
    parse<JSON>::op<Opts>(buff, args...);

    auto to_hictk_chrom = [](const auto& chrom) -> hictk::Chromosome {
      return {chrom.id, chrom.name, chrom.size};
    };

    // clang-format off
    std::ranges::sort(buff, [](const auto& c1, const auto& c2) { return c1.id < c2.id; });
    chroms = buff |
      std::ranges::views::transform(to_hictk_chrom) |
      std::ranges::to<hictk::Reference>();
    // clang-format on
  }
};

template <>
struct to<JSON, hictk::Reference> {
  template <auto Opts>
  static void op(const hictk::Reference& chroms, auto&&... args) noexcept {
    std::vector<const hictk::Chromosome*> buff(chroms.size());
    std::ranges::transform(chroms, buff.begin(), [](const auto& chrom) { return &chrom; });
    serialize<JSON>::op<Opts>(buff, args...);
  }
};

}  // namespace glz

namespace nchg {

template <typename T>
inline std::string to_json_string(const T& value) {
  std::string buff;
  to_json_string(value, buff);

  return buff;
}

template <typename T>
inline void to_json_string(const T& value, std::string& buff) {
  auto ec = glz::write_json(value, buff);
  if (ec) {
    throw std::runtime_error(glz::format_error(ec));
  }
}

template <typename T>
inline T parse_json_string(std::string_view buff) {
  auto res = glz::read_json<T>(buff);
  if (res.has_value()) {
    return res.value();
  }
  throw std::runtime_error(glz::format_error(res.error()));
}

inline glz::json_t parse_json_string(std::string_view buff) {
  return parse_json_string<glz::json_t>(buff);
}

}  // namespace nchg
