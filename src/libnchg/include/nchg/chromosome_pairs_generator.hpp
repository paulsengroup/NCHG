// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>

#pragma once

#ifdef __cpp_lib_generator
#include <generator>
#else
#include <vector>
#endif

#include <hictk/chromosome.hpp>
#include <hictk/reference.hpp>
#include <utility>

namespace nchg {

#ifdef __cpp_lib_generator
[[nodiscard]] inline std::generator<std::pair<hictk::Chromosome, hictk::Chromosome>>
generate_chromosome_pairs_upper_triangle(hictk::Reference reference, bool include_ALL = false) {
  for (const auto& chrom1 : reference) {
    if (chrom1.is_all() && !include_ALL) {
      continue;
    }
    for (auto chrom2_id = chrom1.id(); chrom2_id < reference.size(); ++chrom2_id) {
      const auto& chrom2 = reference.at(chrom2_id);
      if (chrom2.is_all() && !include_ALL) {
        continue;
      }

      co_yield std::make_pair(chrom1, chrom2);
    }
  }
}
#else
[[nodiscard]] constexpr std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>>
generate_chromosome_pairs_upper_triangle(const hictk::Reference& reference,
                                         bool include_ALL = false) {
  std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> pairs{};
  for (const auto& chrom1 : reference) {
    if (chrom1.is_all() && !include_ALL) {
      continue;
    }
    for (auto chrom2_id = chrom1.id(); chrom2_id < reference.size(); ++chrom2_id) {
      const auto& chrom2 = reference.at(chrom2_id);
      if (chrom2.is_all() && !include_ALL) {
        continue;
      }

      pairs.emplace_back(chrom1, chrom2);
    }
  }
  return pairs;
}

#endif

}  // namespace nchg
