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

#pragma once

#include <parallel_hashmap/btree.h>

#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <string>
#include <vector>

#include "nchg/parquet_stats_file_writer.hpp"
#include "nchg/tmpdir.hpp"

namespace nchg {

namespace internal {
// This is defined outside ParquetStatsFileMerger to avoid issues with default initialization of
// nested struct
struct ParquetStatsMergerParams {
  std::filesystem::path tmpdir{TmpDir::default_temp_directory_path()};
  std::uint64_t threads{1};
  std::uint64_t memory_limit_mb{4096};
  std::string compression_method{"zstd"};
  std::uint8_t compression_level{1};
};

}  // namespace internal

class ParquetStatsFileMerger {
  std::vector<std::vector<std::filesystem::path>> _file_groups;

 public:
  using Params = internal::ParquetStatsMergerParams;

  static constexpr auto format_version = ParquetStatsFileWriter::format_version;
  explicit ParquetStatsFileMerger(
      const phmap::btree_map<hictk::Chromosome, std::vector<std::filesystem::path>>& files);

  std::size_t merge(const std::filesystem::path& dest, bool sort, const Params& params = {});
};

}  // namespace nchg
