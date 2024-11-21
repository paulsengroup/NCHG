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

#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <glaze/glaze_exceptions.hpp>
#include <sstream>
#include <string_view>

#include "./common.hpp"
#include "nchg/file_metadata.hpp"
#include "nchg/test/tmpdir.hpp"

namespace nchg::test {

TEST_CASE("FileResultMetadata", "[short][io][file_store]") {
  boost::mt19937_64 rand_eng{3294462971978216056ULL};  // NOLINT

  const auto data_dir = testdir() / "file_result_metadata_001";
  std::filesystem::remove_all(data_dir);        // NOLINT
  std::filesystem::create_directory(data_dir);  // NOLINT

  static const auto cwd = std::filesystem::current_path();
  std::filesystem::current_path(data_dir);

  constexpr std::size_t size1 = 1024;
  constexpr std::size_t size2 = 2048;

  const auto path1 = generate_random_file("file1", rand_eng, size1);
  const auto path2 = generate_random_file("dir/file2", rand_eng, size2);

  const std::string checksum1{"f4e7806963c3626a2c94a450f7c3c8cb"};
  const std::string checksum2{"f3b0311bc7ec3c441ddc557e7de7a984"};

  using FileMetadata = NCHGResultMetadata::FileMetadata;
  const FileMetadata metadata1{path1, XXH3Digest{checksum1}, size1};
  const FileMetadata metadata2{path2, XXH3Digest{checksum2}, size2};

  // clang-format off
  // printf 'NCHG v0.0.2XXH3 (128 bits)5368709121.0dir/file2f3b0311bc7ec3c441ddc557e7de7a9842048file1f4e7806963c3626a2c94a450f7c3c8cb1024' | xxhsum -H128
  constexpr std::string_view expected_digest{"a0563c02213cc9c23c46a75f4097c4c0"};
  // clang-format on

  SECTION("valid") {
    const std::string test_file_content = fmt::format(
        R"~({{
      "format": "NCHG metadata",
      "format-version": "1.0",
      "created-by": "NCHG v0.0.2",
      "creation-time": "2024-11-18T10:55:25",
      "digest": "{}",
      "digest-algorithm": "XXH3 (128 bits)",
      "digest-sample-size": 536870912,
      "records": [
          {{
            "name": "{}",
            "digest": "{}",
            "size": {}
          }},
          {{
            "name": "{}",
            "digest": "{}",
            "size": {}
          }}
      ]
      }})~",
        expected_digest, path1.string(), checksum1, size1, path2.string(), checksum2, size2);

    std::stringstream ss;
    ss << test_file_content;

    SECTION("accessors") {
      const auto metadata = NCHGResultMetadata::from_stream(ss, data_dir);
      CHECK(metadata.path().empty());
      CHECK(metadata.format_version() == "1.0");
      CHECK(metadata.created_by() == "NCHG v0.0.2");
      CHECK(!metadata.creation_time().empty());
      CHECK(metadata.digest() == expected_digest);
      CHECK(metadata.records().size() == 2);

      CHECK(metadata.contains(path1));
      CHECK(metadata.contains(path2));
      CHECK(metadata.at(path1) == metadata1);
      CHECK(metadata.at(path2) == metadata2);

      CHECK_FALSE(metadata.contains(data_dir / "file3"));
      CHECK_THROWS_AS(metadata.at(data_dir / "file3"), std::out_of_range);
    }

    SECTION("checksum") { CHECK(NCHGResultMetadata::checksum(ss) == expected_digest); }

    SECTION("path normalization") {
      const auto abs_path1 = std::filesystem::canonical(path1);
      const auto abs_path2 = std::filesystem::canonical(path2);

      const auto metadata = NCHGResultMetadata::from_stream(ss, data_dir);

      CHECK(metadata.at(std::filesystem::path(".") / path1) == metadata1);
      CHECK(metadata.at(std::filesystem::path(".///") / path1) == metadata1);
      CHECK(metadata.at(abs_path1) == metadata1);
      CHECK(metadata.at(abs_path2) == metadata2);

      std::filesystem::current_path(cwd);
      CHECK(metadata.at(path1) == metadata1);
      CHECK(metadata.at(path2) == metadata2);
    }
  }

  SECTION("invalid") {
    // clang-format off
    const NCHGResultMetadata incomplete_report{
      "",
     "NCHG metadata",
     "1.0",
     "NCHG v0.0.2",
     "2024-11-18T10:55:25",
     "00000000000000000000000000000000",
     "XXH3 (128 bits)",
     512UL << 20UL,
     {metadata1, metadata2}
    };

    const NCHGResultMetadata checksum_mismatch{
      "",
     "NCHG metadata",
     "1.0",
     "NCHG v0.0.2",
     "2024-11-18T10:55:25",
     "9842da6b6a331c68899a524b7fd35e34",
     "XXH3 (128 bits)",
     512UL << 20UL,
     {metadata1, metadata2}
    };

    const NCHGResultMetadata invalid_format{
      "",
      "abc",
      "1.0",
      "NCHG v0.0.2",
      "2024-11-18T10:55:25",
      "a0563c02213cc9c23c46a75f4097c4c0",
      "XXH3 (128 bits)",
      512UL << 20UL,
      {metadata1, metadata2}
    };

    const NCHGResultMetadata invalid_algorithm{
      "",
     "NCHG metadata",
      "1.0",
      "NCHG v0.0.2",
      "2024-11-18T10:55:25",
      "28adcb5a80d965246014d2a20522cbb3",
      "SHA256",
      512UL << 20UL,
      {metadata1, metadata2}
    };

    const NCHGResultMetadata invalid_sample_size{
      "",
      "NCHG metadata",
      "1.0",
      "NCHG v0.0.2",
      "2024-11-18T10:55:25",
      "a0563c02213cc9c23c46a75f4097c4c0",
      "XXH3 (128 bits)",
      0,
      {metadata1, metadata2}
    };

    const NCHGResultMetadata::FileMetadata invalid_file1{
      "not-a-file",
      XXH3Digest{"99aa06d3014798d86001c324468d497f"},
      0
    };

    const NCHGResultMetadata::FileMetadata invalid_file2{
      path1,
      XXH3Digest{checksum2},
      size1
    };

    const NCHGResultMetadata::FileMetadata invalid_file3{
      path1,
      XXH3Digest{checksum1},
      123
    };

    const NCHGResultMetadata invalid_files1{
      "",
      "NCHG metadata",
      "1.0",
      "NCHG v0.0.2",
      "2024-11-18T10:55:25",
      "0d997983eda168c050aaf68582ccc7a9",
      "XXH3 (128 bits)",
      512UL << 20UL,
      {invalid_file1, metadata2}
    };

    const NCHGResultMetadata invalid_files2{
      "",
      "NCHG metadata",
      "1.0",
      "NCHG v0.0.2",
      "2024-11-18T10:55:25",
      "f6c6e84e2f40fcf37973f2d86649697f",
      "XXH3 (128 bits)",
      512UL << 20UL,
      {invalid_file2, metadata2}
    };

    const NCHGResultMetadata invalid_files3{
      "",
      "NCHG metadata",
      "1.0",
      "NCHG v0.0.2",
      "2024-11-18T10:55:25",
      "0f3ced7734232e6f1e95333fcac9670d",
      "XXH3 (128 bits)",
      512UL << 20UL,
      {invalid_file3, metadata2}
    };
    // clang-format on

    std::stringstream ss;

    ss.str(glz::write_json(incomplete_report).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("report was never finalized"));

    ss.str(glz::write_json(checksum_mismatch).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("checksum mismatch"));

    ss.str(glz::write_json(invalid_format).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("unrecognized format"));

    ss.str(glz::write_json(invalid_algorithm).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("unrecognized digest-algorithm"));

    ss.str(glz::write_json(invalid_sample_size).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("digest-sample-size cannot be 0"));

    ss.str(glz::write_json(invalid_files1).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("file does not exist"));

    ss.str(glz::write_json(invalid_files2).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("checksum mismatch"));

    ss.str(glz::write_json(invalid_files3).value_or("error"));
    ss.clear();
    CHECK_THROWS_WITH(NCHGResultMetadata::from_stream(ss, data_dir),
                      Catch::Matchers::ContainsSubstring("file size mismatch"));
  }

  std::filesystem::current_path(cwd);
}

}  // namespace nchg::test
