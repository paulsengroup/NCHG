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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>

namespace nchg::test {

inline std::filesystem::path stage_file(const std::filesystem::path& src,
                                        const std::filesystem::path& dest_dir) {
  assert(std::filesystem::exists(src));
  assert(std::filesystem::is_directory(dest_dir));

  auto dest = dest_dir / src.filename();
  std::filesystem::copy_file(src, dest, std::filesystem::copy_options::overwrite_existing);
  return dest;
}

[[nodiscard]] inline std::filesystem::path generate_random_file(std::filesystem::path dest,
                                                                boost::mt19937_64& rand_eng,
                                                                std::size_t size) {
  assert(size > 0);
  if (dest.has_parent_path()) {
    std::filesystem::create_directories(dest.parent_path());  // NOLINT
  }

  std::vector<std::int8_t> data(size);
  std::ranges::generate(data, [&]() {
    return boost::random::uniform_int_distribution<std::int8_t>{-125, 125}(rand_eng);
  });

  std::ofstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);
  std::ofstream file(dest, std::ios::trunc | std::ios::binary);

  // NOLINTNEXTLINE
  file.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(size));

  return dest;
}

}  // namespace nchg::test
