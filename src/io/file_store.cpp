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

#include "nchg/file_store.hpp"

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/btree.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <hictk/hic.hpp>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "nchg/text.hpp"

namespace nchg {

[[nodiscard]] static std::filesystem::path sanitize_root_path(const std::filesystem::path& path) {
  if (!std::filesystem::is_directory(path)) {
    throw std::runtime_error("root folder does not exist");
  }
  return std::filesystem::absolute(path).make_preferred();
}

[[nodiscard]] static std::filesystem::path validate_report_name(std::filesystem::path path) {
  if (path.has_parent_path()) {
    throw std::runtime_error("report name cannot have a parent path");
  }
  return path;
}

FileStore::FileStore(const std::filesystem::path& folder, bool force,
                     const std::filesystem::path& report_name)
    : _metadata(NCHGResultMetadata::init_empty(sanitize_root_path(folder) /
                                               validate_report_name(report_name))),
      _report_fs(init_report(_metadata, force)) {
  SPDLOG_INFO("initialized file store under prefix \"{}\"", root().string());
}

FileStore::FileStore(FileStore&& other) noexcept
    : _metadata(std::move(other._metadata)),
      _report_fs(std::move(other._report_fs)),
      _finalized(other._finalized) {}

FileStore::~FileStore() noexcept {
  if (finalized()) {
    return;
  }

  try {
    finalize();
    // finalize() catches all exceptions and re-throws them as std::runtime_errors
    // so this should handle all exception cases
  } catch (const std::runtime_error& e) {
    SPDLOG_ERROR(e.what());
    SPDLOG_ERROR("files located under prefix \"{}\" are likely corrupted or incomplete!",
                 root().string());
  }
}

FileStore& FileStore::operator=(FileStore&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  try {
    if (!finalized()) {
      finalize();
    }
    // finalize() catches all exceptions and re-throws them as std::runtime_errors
    // so this should handle all exception cases
  } catch (const std::runtime_error& e) {
    SPDLOG_ERROR(e.what());
    SPDLOG_ERROR("files located under prefix \"{}\" are likely corrupted or incomplete!",
                 root().string());
  }

  _metadata = std::move(other._metadata);
  _report_fs = std::move(other._report_fs);
  _finalized = other._finalized;

  return *this;
}

const std::filesystem::path& FileStore::report_path() const noexcept { return _metadata.path(); }

std::filesystem::path FileStore::root() const { return _metadata.path().parent_path(); }

auto FileStore::at(const std::filesystem::path& path) const -> NCHGResultMetadata::FileMetadata {
  return _metadata.at(path);
}

auto FileStore::get_registered_files() const noexcept
    -> decltype(std::declval<NCHGResultMetadata>().records()) {
  return _metadata.records();
}

bool FileStore::contains(const std::filesystem::path& path) const {
  return _metadata.contains(path);
}

bool FileStore::finalized() const noexcept { return _finalized; }

void FileStore::register_file(const std::filesystem::path& path) {
  auto validate_path = [&] {
    if (contains(path)) {
      throw std::runtime_error(
          fmt::format("file \"{}\" has already been added to the file store!", path));
    }
  };

  SPDLOG_DEBUG("FileStore: registering file \"{}\"...", path.string());

  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error(
        fmt::format("file does not exist or is not a regular file: \"{}\"", path.string()));
  }

  std::unique_lock lck(_mtx);
  validate_path();
  lck.unlock();

  SPDLOG_DEBUG("hashing file \"{}\"...", path.string());
  auto digest = hash_file(path, static_cast<std::streamsize>(_metadata.digest_sample_size()));
  const auto file_size = std::filesystem::file_size(path);
  SPDLOG_DEBUG("{}: XXH3={}", digest, path.string());

  lck.lock();
  validate_path();
  _metadata.add_record(path, XXH3Digest{std::move(digest)}, file_size);
}

static void rename_file(const std::filesystem::path& src, const std::filesystem::path& dest) {
  try {
    std::filesystem::rename(src, dest);
  } catch (const std::filesystem::filesystem_error&) {
    // NOLINTNEXTLINE
    std::filesystem::copy_file(src, dest, std::filesystem::copy_options::overwrite_existing);
    std::filesystem::remove(src);  // NOLINT
  }
}

void FileStore::move_file_and_register(const std::filesystem::path& path,
                                       std::filesystem::path dest) {
  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error(
        fmt::format("file does not exist or is not a regular file: \"{}\"", path.string()));
  }

  if (dest.empty()) {
    dest = root() / path.filename();
  } else {
    assert(!dest.is_absolute());
    dest = root() / dest;
  }

  try {
    if (std::filesystem::exists(dest)) {
      throw std::runtime_error("file already exists");
    }
    if (dest.has_parent_path()) {
      std::filesystem::create_directories(dest.parent_path());  // NOLINT
    }
    rename_file(path, dest);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("failed to move file from \"{}\" to \"{}\": {}", path, dest, e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to move file from \"{}\" to \"{}\": unknown error", path, dest));
  }

  register_file(dest);
}

std::ofstream FileStore::init_report(const NCHGResultMetadata& metadata, bool force) {
  std::ofstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  try {
    SPDLOG_DEBUG("initializing report file at \"{}\"...", metadata.path().string());
    if (force) {
      std::filesystem::remove(metadata.path());  // NOLINT
    } else if (std::filesystem::exists(metadata.path())) {
      throw std::runtime_error("file already exists");
    }

    const auto root_folder = metadata.path().parent_path();

    if (!std::filesystem::exists(root_folder)) {
      std::filesystem::create_directories(root_folder);
    }
#ifdef __cpp_lib_ios_noreplace
    fs.open(metadata.path(), std::ios::out | std::ios::trunc | std::ios::noreplace);
#else
    fs.open(metadata.path(), std::ios::out | std::ios::trunc);
#endif

    std::string buffer{};
    constexpr glz::opts opts{.prettify = true, .indentation_width = 2};
    const auto ec = glz::write<opts>(metadata, buffer);
    if (ec) {
      throw std::runtime_error(glz::format_error(ec));
    }

    fs.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    fs.flush();

    return fs;

  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format("failed to initialize report file at \"{}\": {}",
                                         metadata.path().string(), e.what()));
  } catch (...) {
    throw std::runtime_error(fmt::format(
        "failed to initialize report file at \"{}\": unknown error", metadata.path().string()));
  }
}

void FileStore::finalize() {
  try {
    SPDLOG_DEBUG("finalizing report file at \"{}\"...", report_path().string());
    [[maybe_unused]] const auto lck = std::scoped_lock{_mtx};
    if (_finalized) {
      throw std::runtime_error("FileStore::finalize() has already been called!");
    }

    assert(_report_fs.is_open());
    for (const auto& record : _metadata.records()) {
      auto digest = hash_file(root() / record.name,
                              static_cast<std::streamsize>(_metadata.digest_sample_size()));
      record.validate(digest, root());
    }

    std::string buffer{};
    constexpr glz::opts opts{.prettify = true, .indentation_width = 2};
    const auto ec = glz::write<opts>(_metadata, buffer);
    if (ec) {
      throw std::runtime_error(glz::format_error(ec));
    }

    _report_fs.seekp(0, std::ios::beg);
    _report_fs.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    _report_fs.flush();
    assert(_report_fs.tellp() ==
           static_cast<std::streamsize>(std::filesystem::file_size(report_path())));

    _finalized = true;
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("failed to finalize report file \"{}\": {}", report_path(), e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to finalize report file \"{}\": unknown error", report_path()));
  }
}

}  // namespace nchg
