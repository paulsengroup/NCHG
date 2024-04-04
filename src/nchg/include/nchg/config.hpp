// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#pragma once

#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <variant>

namespace nchg {
struct ComputePvalConfig {
  std::filesystem::path path{};
  std::uint32_t resolution{};
  std::string chrom1{"all"};
  std::string chrom2{"all"};

  bool cis_only{false};
  bool trans_only{false};

  std::uint64_t min_delta{40'000};
  std::uint64_t max_delta{std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t num_quantiles{100};

  bool write_header{true};

  std::uint8_t verbosity{4};
};

struct FilterConfig {
  std::filesystem::path path{};

  double fdr{0.01};
  double log_ratio{2.0};
  bool keep_non_significant{false};

  bool write_header{true};

  std::uint8_t verbosity{4};
};

struct ExpectedConfig {
  std::filesystem::path input_path{};
  std::uint32_t resolution{};

  std::filesystem::path output_path{};
  bool force{false};

  std::string chrom1{"all"};
  std::string chrom2{"all"};

  std::uint64_t min_delta{40'000};
  std::uint64_t max_delta{std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t num_quantiles{100};

  std::uint8_t verbosity{4};
};

using Config = std::variant<std::monostate, ComputePvalConfig, ExpectedConfig, FilterConfig>;

}  // namespace nchg
