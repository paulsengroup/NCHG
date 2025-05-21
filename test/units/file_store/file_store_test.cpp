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

#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <glaze/glaze_exceptions.hpp>
#include <nchg/file_store.hpp>
#include <stdexcept>
#include <string_view>

#include "./common.hpp"
#include "nchg/test/tmpdir.hpp"
#include "nchg/version.hpp"

namespace nchg::test {

static std::filesystem::path stage_file(const std::filesystem::path& src, const FileStore& store) {
  return stage_file(src, store.root());
}

static void validate_report_placeholder(const std::filesystem::path& path) {
  std::string buffer{};
  NCHGResultMetadata data{};
  if (const auto ec = glz::read_file_json(data, path.c_str(), buffer); ec) {
    throw std::runtime_error(glz::format_error(ec, buffer));
  }

  CHECK(data.digest() == "00000000000000000000000000000000");
}

TEST_CASE("FileStore", "[short][io][file_store]") {
  const auto test_file1 = datadir / "ENCFF447ERX.1000000.compartments.bed";
  const auto test_file2 = datadir / "ENCFF447ERX.1000000.cool";

  // xxhsum -H128 myfile
  constexpr std::string_view test_hash1{"5b521558895b0836d2cad90530e970bf"};
  constexpr std::string_view test_hash2{"53723b55cdd873ae74a63493432c7579"};

  // xxhsum -H128 <(cat <(head -c 524288 myfile) <(tail -c 524288 myfile))
  // 1048576 bytes are enough to hash test_file1 in its entirety, but only hash parts of test_file2
  // constexpr std::size_t sample_size_bytes = 1'048'576;
  // constexpr std::string_view test_hash3 = test_hash1;
  // constexpr std::string_view test_hash4{"bd371616784813f444d389668fc14fef"};

  SECTION("Ctor") {
    const auto store_dir = testdir() / "file_store_001";

    CHECK_THROWS_WITH(FileStore(store_dir, false),
                      Catch::Matchers::ContainsSubstring("root folder does not exist"));

    std::filesystem::create_directories(store_dir);  // NOLINT
    CHECK_THROWS_WITH(FileStore(store_dir, false, std::filesystem::path{"parent"} / "report.json"),
                      Catch::Matchers::ContainsSubstring("cannot have a parent path"));

    {
      const FileStore store1(store_dir, false);
      const FileStore store2(store_dir, false, "new_report.json");
      CHECK_THROWS_WITH(FileStore(store_dir, false),
                        Catch::Matchers::ContainsSubstring("file already exists"));
    }
    CHECK_THROWS_WITH(FileStore(store_dir, false),
                      Catch::Matchers::ContainsSubstring("file already exists"));
    CHECK_NOTHROW(FileStore{store_dir, true});
  }

  SECTION("accessors") {
    const auto store_dir = std::filesystem::canonical(testdir()) / "file_store_002";
    std::filesystem::create_directories(store_dir);  // NOLINT

    const FileStore store(store_dir, false);

    CHECK(store.root() == store_dir);
    CHECK(store.report_path() == store_dir / "report.json");
    CHECK(store.get_registered_files().empty());
    CHECK_THROWS_AS(store.at("foo"), std::out_of_range);
    CHECK(!store.contains("foo"));
  }

  SECTION("registration") {
    const auto store_dir = testdir() / "file_store_003";
    std::filesystem::remove_all(store_dir);          // NOLINT
    std::filesystem::create_directories(store_dir);  // NOLINT

    FileStore store(store_dir, false);

    SECTION("invalid paths") {
      CHECK_THROWS_WITH(store.register_file(""),
                        Catch::Matchers::ContainsSubstring("file does not exist"));
      CHECK_THROWS_WITH(store.register_file(store_dir / "not-a-file"),
                        Catch::Matchers::ContainsSubstring("file does not exist"));
      std::filesystem::create_directories(store_dir / "a-folder");  // NOLINT
      CHECK_THROWS_WITH(store.register_file(store_dir / "a-folder"),
                        Catch::Matchers::ContainsSubstring("is not a regular file"));
    }

    SECTION("valid paths") {
      const auto path1 = stage_file(test_file1, store);

      store.register_file(path1);
      CHECK(store.get_registered_files().size() == 1);
      CHECK_THROWS_WITH(store.register_file(store_dir / test_file1.filename()),
                        Catch::Matchers::ContainsSubstring("has already been added"));

      const auto path2 = stage_file(test_file2, store);
      store.register_file(path2);
      CHECK(store.get_registered_files().size() == 2);

      std::filesystem::create_directories(store_dir / "dir");  // NOLINT
      const auto path3 = stage_file(test_file2, store_dir / "dir");
      store.register_file(path3);
      CHECK(store.get_registered_files().size() == 3);

      CHECK(store.at(path1).digest() == test_hash1);
      CHECK(store.at(path2).digest() == test_hash2);
      CHECK(store.at(path3).digest() == test_hash2);

      CHECK_THROWS_AS(store.at("not-a-file"), std::out_of_range);
    }

    SECTION("with move") {
      const auto scratch_dir = testdir() / "scratch";
      std::filesystem::create_directories(scratch_dir);  // NOLINT

      REQUIRE(!std::filesystem::exists(scratch_dir / "not-a-file"));
      CHECK_THROWS_WITH(store.move_file_and_register(scratch_dir / "not-a-file"),
                        Catch::Matchers::ContainsSubstring("file does not exist"));

      std::filesystem::create_directories(scratch_dir / "a-folder");  // NOLINT
      CHECK_THROWS_WITH(store.move_file_and_register(scratch_dir / "a-folder"),
                        Catch::Matchers::ContainsSubstring("is not a regular file"));

      REQUIRE(store.get_registered_files().empty());

      const auto path1 = stage_file(test_file1, scratch_dir);
      store.move_file_and_register(path1);
      CHECK(!std::filesystem::exists(path1));
      CHECK(store.at(store_dir / path1.filename()).digest() == test_hash1);

      const auto path2 = stage_file(test_file1, scratch_dir);
      CHECK_THROWS(store.move_file_and_register(path2),
                   Catch::Matchers::ContainsSubstring("file already exists"));
      CHECK(std::filesystem::exists(path2));

      store.move_file_and_register(path2, "custom-file-name");
      CHECK(!std::filesystem::exists(path2));
      CHECK(store.at(store_dir / "custom-file-name").digest() == test_hash1);

      const auto path3 = stage_file(test_file1, scratch_dir);
      CHECK_THROWS(store.move_file_and_register(path3, "custom-file-name"),
                   Catch::Matchers::ContainsSubstring("has already been added"));
      CHECK(std::filesystem::exists(path3));
    }
  }

  SECTION("report") {
    const auto store_dir = testdir() / "file_store_004";
    std::filesystem::remove_all(store_dir);          // NOLINT
    std::filesystem::create_directories(store_dir);  // NOLINT

    FileStore store(store_dir, false);

    validate_report_placeholder(store.report_path());

    const auto path1 = stage_file(test_file1, store);
    const auto path2 = stage_file(test_file2, store);
    store.register_file(path1);
    store.register_file(path2);

    REQUIRE(store.at(path1).digest() == test_hash1);
    REQUIRE(store.at(path2).digest() == test_hash2);

    store.finalize();
    REQUIRE(store.finalized());

    const auto report = NCHGResultMetadata::from_file(store.report_path());
    // const auto msg = glz::write_json(report);
    // std::ofstream ofs("/tmp/test.json", std::ios::trunc);
    // ofs.write(msg->data(), static_cast<std::streamsize>(msg->size()));

    CHECK(report.path() == store.report_path());
    CHECK(report.records().size() == 2);
    CHECK(report.at(path1.filename()).digest() == test_hash1);
    CHECK(report.at(path2.filename()).digest() == test_hash2);

    CHECK(report.format_version() == "1.0");
    CHECK(report.created_by() == config::version::str_long());
    CHECK(report.digest().size() == 32);
    CHECK(report.digest_algorithm() == "XXH3 (128 bits)");
    CHECK(report.digest_sample_size() == 512UL << 20UL);
  }
}

}  // namespace nchg::test
