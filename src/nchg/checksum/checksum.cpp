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

#include <spdlog/spdlog.h>

#include <chrono>

#include "nchg/file_hashing.hpp"
#include "nchg/file_metadata.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

int run_command(const ChecksumConfig& c) {
  const auto t0 = std::chrono::steady_clock::now();
  for (const auto& p : c.files) {
    if (p.extension() == ".json") {
      try {
        const auto metadata = NCHGResultMetadata::from_file(p, false);
        const auto digest = metadata.checksum();
        if (metadata.digest() != digest) {
          SPDLOG_WARN(
              "checksum for \"{}\" did not match the precomputed checksum: expected {}, found {}",
              p.string(), metadata.digest(), digest);
        }
        fmt::print("{}  {}\n", p.string(), digest);
        continue;
        // NOLINTNEXTLINE
      } catch (...) {
      }
    }
    const auto digest = hash_file(p, 512UL << 20UL);
    fmt::print("{}  {}\n", p.string(), digest);
  }

  const auto t1 = std::chrono::steady_clock::now();
  if (c.files.size() == 1) {
    SPDLOG_INFO("checksummed one file in {}", format_duration(t1 - t0));
  } else {
    SPDLOG_INFO("checksummed {} files in {}", c.files.size(), format_duration(t1 - t0));
  }

  return 0;
}

}  // namespace nchg
