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

#include "nchg/file_hashing.hpp"

#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <string_view>

#include "nchg/test/tmpdir.hpp"

namespace nchg::test {

TEST_CASE("File hashing", "[short][io][file_hashing]") {
  const auto test_file1 = datadir / "ENCFF447ERX.1000000.compartments.bed";
  const auto test_file2 = datadir / "ENCFF447ERX.minified.mcool";

  SECTION("full") {
    // xxhsum -H128 myfile
    constexpr std::string_view test_hash1{"5b521558895b0836d2cad90530e970bf"};
    constexpr std::string_view test_hash2{"154545297622748d3dc1608a11b256c1"};

    const auto test_file1_size =
        static_cast<std::streamsize>(std::filesystem::file_size(test_file1));
    const auto test_file2_size =
        static_cast<std::streamsize>(std::filesystem::file_size(test_file2));

    CHECK(hash_file(test_file1, test_file1_size) == test_hash1);
    CHECK(hash_file(test_file2, test_file2_size) == test_hash2);
  }

  SECTION("partial") {
    // 1048576 bytes are enough to hash test_file1 in its entirety, but only hash parts of
    // test_file2
    constexpr std::streamsize sample_size_bytes{1'048'576};

    // xxhsum -H128 myfile
    const std::string_view test_hash1{"5b521558895b0836d2cad90530e970bf"};
    // xxhsum -H128 <(cat <(head -c 524288 myfile) <(tail -c 524288 myfile))
    const std::string_view test_hash2{"b571ab50cd0c1ed7d2ba9b98930d12d7"};

    CHECK(hash_file(test_file1, sample_size_bytes) == test_hash1);
    CHECK(hash_file(test_file2, sample_size_bytes) == test_hash2);
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace nchg::test
