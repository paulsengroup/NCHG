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
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/btree.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <BS_thread_pool.hpp>
#include <compare>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <glaze/glaze_exceptions.hpp>
#include <istream>
#include <string>
#include <string_view>
#include <vector>

#include "nchg/file_hashing.hpp"

namespace nchg {

class NCHGResultMetadata {
 public:
  struct FileMetadata {
    std::filesystem::path name{};
    XXH3Digest digest{};
    std::size_t size{};

    bool operator==(const std::filesystem::path& other) const noexcept;
    bool operator==(const FileMetadata& other) const noexcept;

    std::strong_ordering operator<=>(const FileMetadata& other) const noexcept = default;
    std::weak_ordering operator<=>(const std::filesystem::path& other) const noexcept;

    void validate(std::string_view expected_digest,
                  const std::filesystem::path& root_dir = "") const;

    struct glaze {
      friend struct FileMetadata;
      using T = FileMetadata;
      static constexpr auto value = glz::object(&T::name, &T::digest, &T::size);
    };
  };

  struct ValidationResult {
    bool successfully_finalized{};
    std::string expected_checksum{};
    std::string computed_checksum{};
    std::vector<std::string> report_validation_failures{};
    std::vector<std::pair<std::filesystem::path, std::string>> record_validation_failures{};
    std::size_t num_records{};
    std::exception_ptr unhandled_exception{};

    explicit operator bool() const noexcept;
    [[noreturn]] void throw_exception() const;
  };

 private:
  struct FileMetadataCmp {
    using is_transparent = void;
    bool operator()(const FileMetadata& lhs, const FileMetadata& rhs) const noexcept;
    bool operator()(const FileMetadata& lhs, const std::filesystem::path& rhs) const noexcept;
    bool operator()(const std::filesystem::path& lhs, const FileMetadata& rhs) const noexcept;
  };

  std::filesystem::path _root_dir{};
  std::filesystem::path _path{};
  std::string _format{"NCHG metadata"};
  std::string _format_version{"1.0"};
  std::string _created_by{"NCHG v0.0.2"};  // TODO fixme
  std::string _creation_time{};
  XXH3Digest _digest{};
  std::string _digest_algorithm{"XXH3 (128 bits)"};
  std::size_t _digest_sample_size{512UL << 20UL};
  phmap::btree_set<FileMetadata, FileMetadataCmp> _records{};

 public:
  NCHGResultMetadata() = default;
  // Note that this ctor does not do any validation on its params
  NCHGResultMetadata(const std::filesystem::path& path_, std::string format_,
                     std::string format_version_, std::string created_by_,
                     std::string creation_time_, std::string digest_, std::string digest_algorithm_,
                     std::size_t digest_sample_size_, const std::vector<FileMetadata>& records_);
  static NCHGResultMetadata init_empty(const std::filesystem::path& path_);
  static NCHGResultMetadata from_file(const std::filesystem::path& path_, bool validate_ = true);
  static NCHGResultMetadata from_stream(std::istream& stream, const std::filesystem::path& root_dir,
                                        bool validate_ = true);

  [[nodiscard]] auto validate(BS::thread_pool* tpool = nullptr) const noexcept -> ValidationResult;
  void validate_or_throw() const;
  [[nodiscard]] std::string checksum() const;
  [[nodiscard]] static std::string checksum(const std::filesystem::path& path_);
  [[nodiscard]] static std::string checksum(std::istream& stream);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string_view format_version() const noexcept;
  [[nodiscard]] std::string_view created_by() const noexcept;
  [[nodiscard]] std::string_view creation_time() const noexcept;
  [[nodiscard]] std::string_view digest() const noexcept;
  [[nodiscard]] std::string_view digest_algorithm() const noexcept;
  [[nodiscard]] std::size_t digest_sample_size() const noexcept;
  [[nodiscard]] auto records() const noexcept
      -> const phmap::btree_set<FileMetadata, FileMetadataCmp>&;
  [[nodiscard]] auto at(const std::filesystem::path& name) const -> const FileMetadata&;
  [[nodiscard]] bool contains(const std::filesystem::path& name) const;

  void add_record(std::filesystem::path path, XXH3Digest digest, std::size_t size);

  struct glaze {
    friend class NCHGResultMetadata;
    using T = NCHGResultMetadata;
    // clang-format off
    static constexpr auto value =
        glz::object(
          "format", &T::_format,
          "format-version", &T::_format_version,
          "created-by", &T::_created_by,
          "creation-time", &T::_creation_time,
          "digest", &T::_digest,
          "digest-algorithm", &T::_digest_algorithm,
          "digest-sample-size", &T::_digest_sample_size,
          "records", &T::_records
        );
    // clang-format on
  };

 private:
  [[nodiscard]] static XXH3Digest null_digest();
  [[nodiscard]] static std::string current_time();
};

}  // namespace nchg
