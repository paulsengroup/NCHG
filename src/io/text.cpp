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

#include "nchg/text.hpp"

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <hictk/chromosome.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <istream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nchg {

[[nodiscard]] static phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, std::istream &stream,
    const std::filesystem::path &path) {
  if (path.empty()) {
    SPDLOG_INFO("reading the user-provided bin mask...");
  } else {
    SPDLOG_INFO("reading the user-provided bin mask from {}...", path);
  }

  phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> mask{};
  std::string buffer{};

  std::size_t i = 1;
  try {
    for (; std::getline(stream, buffer); ++i) {
      if (buffer.empty()) {
        continue;
      }

      if (buffer.back() == '\r') {
        buffer.resize(buffer.size() - 1);
      }

      try {
        const auto record = truncate_record<3>(buffer);
        const auto domain = hictk::GenomicInterval::parse_bed(chroms, record);

        const auto num_bins = (domain.chrom().size() + bin_size - 1) / bin_size;
        auto [match, _] = mask.try_emplace(domain.chrom(), std::vector(num_bins, false));

        const std::size_t j0 = domain.start() / bin_size;
        const std::size_t j1 = (domain.end() / bin_size) + 1;

        for (std::size_t j = j0; j < j1; ++j) {
          match->second[j] = true;
        }

      } catch (const std::exception &e) {
        if (path.empty()) {
          throw std::runtime_error(fmt::format(
              "found an invalid record at line {} while parsing the bon mask: {}", i, e.what()));
        }
        throw std::runtime_error(
            fmt::format("found an invalid record at line {} of file {}: {}", i, path, e.what()));
      }
    }

  } catch (const std::exception &) {
    if (!stream.eof()) {
      throw;
    }
  }

  std::size_t num_bad_bins = 0;
  for (const auto &v : std::ranges::views::values(mask)) {
    num_bad_bins += std::accumulate(v.begin(), v.end(), 0UZ);
  }

  SPDLOG_INFO("masked {} bad bins based on {} domains read from {}...", num_bad_bins, i - 1, path);
  return mask;
}

phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, const std::string &payload) {
  if (payload.empty()) {
    return {};
  }

  std::stringstream stream(payload);
  return parse_bin_mask(chroms, bin_size, stream, {});
}

phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> parse_bin_mask(
    const hictk::Reference &chroms, std::uint32_t bin_size, const std::filesystem::path &path) {
  if (path.empty()) {
    return {};
  }

  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);
  fs.open(path);
  return parse_bin_mask(chroms, bin_size, fs, path);
}

}  // namespace nchg
