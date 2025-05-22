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

#include <filesystem>
#include <fstream>
#include <mutex>

#include "nchg/file_metadata.hpp"

namespace nchg {

class FileStore {
  NCHGResultMetadata _metadata;
  mutable std::mutex _mtx;
  std::ofstream _report_fs;
  bool _finalized{false};

 public:
  FileStore() = delete;
  FileStore(const std::filesystem::path& folder, bool force,
            const std::filesystem::path& report_name = "report.json");
  FileStore(const FileStore&) = delete;
  FileStore(FileStore&& other) noexcept;

  ~FileStore() noexcept;

  FileStore& operator=(const FileStore&) = delete;
  FileStore& operator=(FileStore&& other) noexcept;

  [[nodiscard]] const std::filesystem::path& report_path() const noexcept;
  [[nodiscard]] std::filesystem::path root() const;

  [[nodiscard]] auto at(const std::filesystem::path& path) const
      -> NCHGResultMetadata::FileMetadata;
  [[nodiscard]] auto get_registered_files() const noexcept
      -> decltype(std::declval<NCHGResultMetadata>().records());
  [[nodiscard]] bool contains(const std::filesystem::path& path) const;
  [[nodiscard]] bool finalized() const noexcept;

  // Paths to be registered must point to regular files and reside under the root path
  void register_file(const std::filesystem::path& path);
  // File pointed by path will be moved under root.
  // When dest is not provided, the file will be placed directly under root with the same filename
  // component as what is found in path
  void move_file_and_register(const std::filesystem::path& path, std::filesystem::path dest = "");

  // It is usually not necessary to manually call finalize, as by default the file store is
  // finalized when the destructor is called
  void finalize();

 private:
  [[nodiscard]] static std::ofstream init_report(const NCHGResultMetadata& metadata, bool force);
};

}  // namespace nchg
