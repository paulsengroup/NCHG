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

#include "nchg/metadata.hpp"

#include <arrow/io/api.h>
#include <fmt/format.h>
#include <parquet/arrow/reader.h>

#include <exception>
#include <filesystem>
#include <glaze/json.hpp>
#include <stdexcept>
#include <string>

#include "nchg/parquet_stats_file_reader.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

[[nodiscard]] static std::string fetch_raw_metadata(const std::filesystem::path& path) {
  return ParquetStatsFileReader::read_metadata(path);
}

[[nodiscard]] static std::string fetch_metadata(const std::filesystem::path& path) {
  auto metadata = parse_json_string(ParquetStatsFileReader::read_metadata(path));
  std::string buff;

  const auto ec = glz::write<glz::opts{.prettify = true}>(metadata, buff);
  if (ec) {
    throw std::runtime_error(glz::format_error(ec));
  }

  return buff;
}

static void validate_parquet_file(const std::filesystem::path& path) {
  try {
    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(path))
    const auto result = parquet::arrow::OpenFile(fp, arrow::default_memory_pool());
    if (!result.ok()) {
      throw std::runtime_error(result.status().ToString());
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("failed to open file \"{}\" for reading: {}", path.string(), e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to open file \"{}\" for reading: unknown error", path.string()));
  }
}

int run_command(const MetadataConfig& c) {
  try {
    validate_parquet_file(c.input_path);

    const auto metadata = c.raw ? fetch_raw_metadata(c.input_path) : fetch_metadata(c.input_path);

    fmt::println("{}", metadata);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("failed to read metadata from \"{}\": {}", c.input_path.string(), e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to read metadata from \"{}\": unknown error", c.input_path.string()));
  }

  return 0;
}

}  // namespace nchg
