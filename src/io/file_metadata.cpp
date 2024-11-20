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

#include "nchg/file_metadata.hpp"

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/btree.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <glaze/glaze_exceptions.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "nchg/file_hashing.hpp"
#include "nchg/text.hpp"

namespace nchg {

[[nodiscard]] static std::filesystem::path normalize_path(
    const std::filesystem::path& path, const std::filesystem::path& root_dir = {}) {
  auto normalized_path = std::filesystem::weakly_canonical(path).make_preferred();
  if (!root_dir.empty()) {
    normalized_path = std::filesystem::relative(normalized_path, root_dir);
  }

  if constexpr (std::filesystem::path::preferred_separator == '/') {
    return normalized_path;
  }

  auto str = normalized_path.string();
  std::ranges::replace(str, std::filesystem::path::preferred_separator, '/');

  return {str};
}

[[nodiscard]] static std::ifstream open_text_file_checked(const std::filesystem::path& path) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("file does not exist");
  }

  std::ifstream ifs{};
  ifs.exceptions(ifs.exceptions() | std::ios::badbit | std::ios::failbit);
  ifs.open(path);
  return ifs;
}

bool NCHGResultMetadata::FileMetadata::operator==(
    const std::filesystem::path& other) const noexcept {
  return name == other;
}

bool NCHGResultMetadata::FileMetadata::operator==(const FileMetadata& other) const noexcept {
  return name == other.name && digest == other.digest && size == other.size;
}

std::weak_ordering NCHGResultMetadata::FileMetadata::operator<=>(
    const std::filesystem::path& other) const noexcept {
  return name <=> other;
}

void NCHGResultMetadata::FileMetadata::validate(std::string_view expected_digest,
                                                const std::filesystem::path& root_dir) const {
  if (name.empty()) {
    throw std::runtime_error("name cannot be empty");
  }

  const auto path = root_dir.empty() ? name : root_dir / name;
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("file does not exist");
  }

  const auto file_size = std::filesystem::file_size(path);
  if (file_size != size) {
    throw std::runtime_error(
        fmt::format("file size mismatch: expected {}, found {}", size, file_size));
  }

  if (!expected_digest.empty() && expected_digest != digest()) {
    throw std::runtime_error(
        fmt::format("checksum mismatch: expected {}, found {}", digest(), expected_digest));
  }
}

NCHGResultMetadata::NCHGResultMetadata(const std::filesystem::path& path_, std::string format_,
                                       std::string format_version_, std::string created_by_,
                                       std::string creation_time_, std::string digest_,
                                       std::string digest_algorithm_,
                                       std::size_t digest_sample_size_,
                                       const std::vector<FileMetadata>& records_)
    : _root_dir(normalize_path(path_.parent_path())),
      _path(std::filesystem::weakly_canonical(path_)),
      _format(std::move(format_)),
      _format_version(std::move(format_version_)),
      _created_by(std::move(created_by_)),
      _creation_time(std::move(creation_time_)),
      _digest(XXH3Digest{std::move(digest_)}),
      _digest_algorithm(std::move(digest_algorithm_)),
      _digest_sample_size(digest_sample_size_),
      _records(records_.begin(), records_.end()) {}

NCHGResultMetadata NCHGResultMetadata::init_empty(const std::filesystem::path& path_) {
  NCHGResultMetadata metadata{};

  metadata._root_dir = normalize_path(path_).parent_path();
  metadata._path = normalize_path(path_);
  metadata._creation_time = current_time();
  metadata._digest = null_digest();

  return metadata;
}

NCHGResultMetadata NCHGResultMetadata::from_file(const std::filesystem::path& path_) {
  try {
    auto ifs = open_text_file_checked(path_);
    auto metadata = from_stream(ifs, path_.parent_path());
    metadata._path = normalize_path(path_);
    metadata._root_dir = metadata._path.parent_path();
    return metadata;
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("failed to import metadata from \"{}\": {}", path_.string(), e.what()));
  }
}

NCHGResultMetadata NCHGResultMetadata::from_stream(std::istream& stream,
                                                   const std::filesystem::path& root_dir) {
  stream.seekg(0, std::ios::end);
  std::string strbuf{};
  strbuf.resize(static_cast<std::size_t>(stream.tellg()));
  stream.seekg(0, std::ios::beg);
  stream.read(strbuf.data(), static_cast<std::streamoff>(strbuf.size()));

  constexpr glz::opts opts{.error_on_missing_keys = true};

  NCHGResultMetadata data{};
  if (const auto ec = glz::read<opts>(data, strbuf); ec) {
    throw std::runtime_error(glz::format_error(ec, strbuf));
  }

  data._root_dir = root_dir;
  data.validate();
  return data;
}

bool NCHGResultMetadata::FileMetadataCmp::operator()(const FileMetadata& lhs,
                                                     const FileMetadata& rhs) const noexcept {
  return lhs.name < rhs.name;
}

bool NCHGResultMetadata::FileMetadataCmp::operator()(
    const FileMetadata& lhs, const std::filesystem::path& rhs) const noexcept {
  return lhs.name < rhs;
}

bool NCHGResultMetadata::FileMetadataCmp::operator()(const std::filesystem::path& lhs,
                                                     const FileMetadata& rhs) const noexcept {
  return lhs < rhs.name;
}

void NCHGResultMetadata::validate() const {
  if (_format != "NCHG metadata") {
    throw std::runtime_error(fmt::format("unrecognized format \"{}\"", _format));
  }

  if (_digest_algorithm != "XXH3 (128 bits)") {
    throw std::runtime_error(
        fmt::format("unrecognized digest-algorithm \"{}\"", _digest_algorithm));
  }

  if (_digest_sample_size == 0) {
    throw std::runtime_error("digest-sample-size cannot be 0");
  }

  auto computed_digest = checksum();
  if (computed_digest != _digest()) {
    throw std::runtime_error(
        fmt::format("checksum mismatch: expected {}, found {}", _digest(), computed_digest));
  }

  std::size_t i = 0;
  auto first = _records.begin();
  const auto last = _records.end();

  if (first == last) {
    return;
  }

  try {
    while (first != last) {
      const auto path = _root_dir / first->name;
      if (!std::filesystem::exists(path)) {
        throw std::runtime_error("file does not exist");
      }
      computed_digest = hash_file(path, static_cast<std::streamsize>(_digest_sample_size));
      first->validate(computed_digest, _root_dir);
      ++i;
      ++first;
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format("failed to validate record #{} (file \"{}\"): {}", i + 1,
                                         first->name.string(), e.what()));
  }
}

std::string NCHGResultMetadata::checksum() const {
  auto state = init_xxh3_state();
  // Fields should be hashed in alphabetic order, except for the records field, which should
  // always come last
  update_state(state, _created_by);
  update_state(state, _digest_algorithm);
  update_state(state, fmt::to_string(_digest_sample_size));
  update_state(state, _format_version);

  for (const auto& record : _records) {
    update_state(state, record.name.string());
    update_state(state, record.digest());
    update_state(state, fmt::to_string(record.size));
  }

  return to_hex(state);
}

std::string NCHGResultMetadata::checksum(const std::filesystem::path& path_) {
  std::ifstream ifs{};
  ifs.exceptions(ifs.exceptions() | std::ios::badbit | std::ios::failbit);

  try {
    ifs.open(path_);
    return checksum(ifs);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format("failed to hash report file \"{}\": {}", path_, e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to hash report file \"{}\": unknown error", path_));
  }
}

std::string NCHGResultMetadata::checksum(std::istream& stream) {
  auto state = init_xxh3_state();
  stream.seekg(0, std::ios::end);
  std::string strbuf{};
  strbuf.resize(static_cast<std::size_t>(stream.tellg()));
  stream.seekg(0, std::ios::beg);
  stream.read(strbuf.data(), static_cast<std::streamoff>(strbuf.size()));

  constexpr glz::opts opts{.error_on_missing_keys = true};
  NCHGResultMetadata data{};
  if (const auto ec = glz::read<opts>(data, strbuf); ec) {
    throw std::runtime_error(glz::format_error(ec, strbuf));
  }

  data.validate();
  return data.checksum();
}

const std::filesystem::path& NCHGResultMetadata::path() const noexcept { return _path; }

std::string_view NCHGResultMetadata::format_version() const noexcept { return _format_version; }

std::string_view NCHGResultMetadata::created_by() const noexcept { return _created_by; }

std::string_view NCHGResultMetadata::creation_time() const noexcept { return _creation_time; }

std::string_view NCHGResultMetadata::digest() const noexcept { return _digest(); }

std::string_view NCHGResultMetadata::digest_algorithm() const noexcept { return _digest_algorithm; }
std::size_t NCHGResultMetadata::digest_sample_size() const noexcept { return _digest_sample_size; }

auto NCHGResultMetadata::records() const noexcept
    -> const phmap::btree_set<FileMetadata, FileMetadataCmp>& {
  return _records;
}

auto NCHGResultMetadata::at(const std::filesystem::path& name) const -> const FileMetadata& {
  auto it = _records.find(name);
  if (it != _records.end()) {
    return *it;
  }
  it = _records.find(normalize_path(name, _root_dir));
  if (it != _records.end()) {
    return *it;
  }

  throw std::out_of_range(fmt::format("no records found for \"{}\"", name));
}

bool NCHGResultMetadata::contains(const std::filesystem::path& name) const {
  auto it = _records.find(name);
  if (it != _records.end()) {
    return true;
  }

  it = _records.find(normalize_path(name, _root_dir));
  return it != _records.end();
}

void NCHGResultMetadata::add_record(std::filesystem::path path, XXH3Digest digest,
                                    std::size_t size) {
  FileMetadata record{normalize_path(std::move(path), _root_dir), std::move(digest), size};
  record.validate("", _root_dir);
  if (_records.contains(record.name)) {
    throw std::runtime_error(fmt::format("entry for \"{}\" already exists", record.name));
  }
  _records.emplace(std::move(record));
  _digest = XXH3Digest{checksum()};
}

XXH3Digest NCHGResultMetadata::null_digest() {
  return XXH3Digest{"00000000000000000000000000000000"};
}

std::string NCHGResultMetadata::current_time() {
  return fmt::format("{:%FT%T}", fmt::gmtime(std::time(nullptr)));
}

}  // namespace nchg
