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

#include <algorithm>
#include <cstdint>
#include <functional>
#include <hictk/chromosome.hpp>
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

  // NOLINTBEGIN(*-const-correctness)
  std::vector margs1(num_bins1, 0.0);
  std::vector margs2(num_bins2, 0.0);
  // NOLINTEND(*-const-correctness)

  for (const auto& p : pixels) {
    margs1[p.coords.bin1.rel_id()] += conditional_static_cast<double>(p.count);
    margs2[p.coords.bin2.rel_id()] += conditional_static_cast<double>(p.count);

    if (chrom1 == chrom2 && p.coords.bin1.id() != p.coords.bin2.id()) {
      margs1[p.coords.bin2.rel_id()] += conditional_static_cast<double>(p.count);
      margs2[p.coords.bin1.rel_id()] += conditional_static_cast<double>(p.count);
    }
  }

  return {std::move(margs1), std::move(margs2)};
}

[[nodiscard]] inline std::size_t count_bad_bins(const std::vector<bool>& mask) noexcept {
  return static_cast<std::size_t>(std::ranges::count_if(mask, [](const auto x) { return !!x; }));
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
    return {std::vector(num_bins1, false), std::vector(num_bins2, false)};
  }
  auto [margs1, margs2] = internal::compute_marginals(pixels, chrom1, chrom2, resolution);

  auto mask1 = mad_max_filtering(margs1, mad_max);
  auto mask2 = chrom1 == chrom2 ? mask1 : mad_max_filtering(margs2, mad_max);

  [[maybe_unused]] const auto bad_bins1 = internal::count_bad_bins(mask1);
  [[maybe_unused]] const auto bad_bins2 = internal::count_bad_bins(mask2);

  if (bad_bins1 == mask1.size() && bad_bins2 == mask2.size()) {
    SPDLOG_INFO("[{}:{}]: MAD-max masking procedure flagged the entire matrix", chrom1.name(),
                chrom2.name());
  } else {
    SPDLOG_INFO("[{}:{}]: MAD-max masking procedure flagged {}/{} bins for {} and {}/{} for {}",
                chrom1.name(), chrom2.name(), bad_bins1, mask1.size(), chrom1.name(), bad_bins2,
                mask2.size(), chrom2.name());
  }

  return {std::move(mask1), std::move(mask2)};
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline std::vector<bool> mad_max_filtering(const Pixels& pixels, const hictk::Chromosome& chrom,
                                           std::uint32_t resolution, double mad_max) {
  if (mad_max == 0) {
    const auto num_bins = (chrom.size() + resolution - 1) / resolution;
    return std::vector<bool>(num_bins, false);  // NOLINT
  }
  auto margs = internal::compute_marginals(pixels, chrom, chrom, resolution).first;

  auto mask = mad_max_filtering(margs, mad_max);
  [[maybe_unused]] const auto bad_bins = internal::count_bad_bins(mask);
  if (mask.size() == bad_bins) {
    SPDLOG_INFO("[{}]: MAD-max masking procedure flagged the entire matrix", chrom.name());
  } else {
    SPDLOG_INFO("[{}]: MAD-max masking procedure flagged {}/{} bins", chrom.name(), bad_bins,
                mask.size());
  }
  return mask;
}
}  // namespace nchg
