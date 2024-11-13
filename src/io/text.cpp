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

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include "nchg/text.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <hictk/chromosome.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nchg {

phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, const std::filesystem::path &path) {
  if (path.empty()) {
    return {};
  }

  SPDLOG_INFO("reading the user-provided bin mask from {}...", path);
  phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> mask{};
  std::string buffer{};

  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  std::size_t i = 1;
  try {
    fs.open(path);

    for (; std::getline(fs, buffer); ++i) {
      if (buffer.empty()) {
        continue;
      }

      if (buffer.back() == '\r') {
        buffer.resize(buffer.size() - 1);
      }

      try {
        const auto record = truncate_record<3>(buffer);
        auto domain = hictk::GenomicInterval::parse_bed(chroms, record);

        const auto num_bins = (domain.chrom().size() + bin_size - 1) / bin_size;
        auto match = mask.try_emplace(domain.chrom(), std::vector<bool>(num_bins, false));

        const std::size_t j0 = domain.start() / bin_size;
        const std::size_t j1 = (domain.end() / bin_size) + 1;

        for (std::size_t j = j0; j < j1; ++j) {
          match.first->second[j] = true;
        }

      } catch (const std::exception &e) {
        throw std::runtime_error(
            fmt::format("found an invalid record at line {} of file {}: {}", i, path, e.what()));
      }
    }

  } catch (const std::exception &) {
    if (!fs.eof()) {
      throw;
    }
  }

  std::size_t num_bad_bins = 0;
  for (const auto &[_, v] : mask) {
    num_bad_bins += std::accumulate(v.begin(), v.end(), 0uz);
  }

  SPDLOG_INFO("masked {} bad bins based on {} domains read from {}...", num_bad_bins, i - 1, path);
  return mask;
}

}  // namespace nchg
