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

#include <spdlog/spdlog.h>

#include <cstdint>
#include <functional>
#include <hictk/chromosome.hpp>
#include <numeric>
#include <utility>
#include <vector>

#include "nchg/median.hpp"

namespace nchg {
namespace internal {

template <typename Pixels>
  requires PixelRange<Pixels>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>> compute_marginals(
    const Pixels& pixels, const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
    std::uint32_t resolution) {
  const auto num_bins1 = (chrom1.size() + resolution - 1) / resolution;
  const auto num_bins2 = (chrom2.size() + resolution - 1) / resolution;

  std::vector<double> margs1(num_bins1, 0);
  std::vector<double> margs2(num_bins2, 0);

  for (const auto& p : pixels) {
    margs1[p.coords.bin1.rel_id()] += conditional_static_cast<double>(p.count);
    margs2[p.coords.bin2.rel_id()] += conditional_static_cast<double>(p.count);

    if (chrom1 == chrom2 && p.coords.bin1.id() != p.coords.bin2.id()) {
      margs1[p.coords.bin2.rel_id()] += conditional_static_cast<double>(p.count);
      margs2[p.coords.bin1.rel_id()] += conditional_static_cast<double>(p.count);
    }
  }

  return std::make_pair(margs1, margs2);
}
}  // namespace internal

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::pair<std::vector<bool>, std::vector<bool>> mad_max_filtering(
    const Pixels& pixels, const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
    std::uint32_t resolution, double mad_max) {
  if (mad_max == 0) {
    const auto num_bins1 = (chrom1.size() + resolution - 1) / resolution;
    const auto num_bins2 = (chrom2.size() + resolution - 1) / resolution;
    return std::make_pair(std::vector<bool>(num_bins1, false), std::vector<bool>(num_bins2, false));
  }
  auto [margs1, margs2] = internal::compute_marginals(pixels, chrom1, chrom2, resolution);

  auto mask1 = mad_max_filtering(margs1, mad_max);
  auto mask2 = chrom1 == chrom2 ? mask1 : mad_max_filtering(margs2, mad_max);

  SPDLOG_INFO(
      FMT_STRING("[{}:{}]: MAD-max masking procedure flagged {}/{} bins for {} and {}/{} for {}"),
      chrom1.name(), chrom2.name(), std::ranges::fold_left(mask1, 0, std::plus{}), mask1.size(),
      chrom1.name(), std::accumulate(mask2.begin(), mask2.end(), 0), mask2.size(), chrom2.name());

  return std::make_pair(mask1, mask2);
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::vector<bool> mad_max_filtering(const Pixels& pixels, const hictk::Chromosome& chrom,
                                           std::uint32_t resolution, double mad_max) {
  if (mad_max == 0) {
    const auto num_bins = (chrom.size() + resolution - 1) / resolution;
    return std::vector<bool>(num_bins, false);  // NOLINT
  }
  auto [margs, _] = internal::compute_marginals(pixels, chrom, chrom, resolution);

  auto mask = mad_max_filtering(margs, mad_max);
  SPDLOG_INFO(FMT_STRING("[{}]: MAD-max masking procedure flagged {}/{} bins"), chrom.name(),
              std::accumulate(mask.begin(), mask.end(), 0), mask.size());
  return mask;
}
}  // namespace nchg
