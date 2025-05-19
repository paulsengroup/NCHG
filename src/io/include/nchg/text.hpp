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

#include <parallel_hashmap/phmap.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/reference.hpp>
#include <string_view>
#include <vector>

namespace nchg {

template <std::size_t NTOKS>
[[nodiscard]] constexpr std::string_view truncate_record(std::string_view record, char sep = '\t');

template <std::size_t NTOKS>
[[nodiscard]] constexpr std::array<std::string_view, NTOKS> tokenize_record(std::string_view record,
                                                                            char sep = '\t');

[[nodiscard]] phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, const std::filesystem::path &path);

}  // namespace nchg

#include "../../text_impl.hpp"
