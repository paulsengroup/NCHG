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
#include <cassert>
#include <cmath>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>
#include <stdexcept>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/concepts.hpp"

namespace nchg {
namespace internal {
template <typename T>
[[nodiscard]] inline T median(std::vector<T> v) {
  if (v.empty()) {
    throw std::runtime_error("median was called on an empty vector");
  }

  const auto size = static_cast<std::ptrdiff_t>(v.size());
  auto first = v.begin();
  auto mid = first + (size / 2);
  auto last = v.end();

  std::nth_element(first, mid, last);

  if (size % 2 != 0) {
    return *mid;
  }

  const auto n1 = *mid;
  std::nth_element(first, --mid, last);
  const auto n2 = *mid;

  return (n1 + n2) / 2;
}

template <typename PixelIt>
[[nodiscard]] inline std::pair<std::vector<double>, std::vector<double>> compute_marginals(
    PixelIt first_pixel, PixelIt last_pixel, const hictk::Chromosome& chrom1,
    const hictk::Chromosome& chrom2, std::uint32_t resolution) {
  using N = decltype(first_pixel->count);

  const auto num_bins1 = (chrom1.size() + resolution - 1) / resolution;
  const auto num_bins2 = (chrom2.size() + resolution - 1) / resolution;

  std::vector<double> margs1(num_bins1, 0);
  std::vector<double> margs2(num_bins2, 0);

  std::for_each(first_pixel, last_pixel, [&](const hictk::Pixel<N>& p) {
    margs1[p.coords.bin1.rel_id()] += conditional_static_cast<N>(p.count);
    margs2[p.coords.bin2.rel_id()] += conditional_static_cast<N>(p.count);

    if (chrom1 == chrom2 && p.coords.bin1.id() != p.coords.bin2.id()) {
      margs1[p.coords.bin2.rel_id()] += conditional_static_cast<N>(p.count);
      margs2[p.coords.bin1.rel_id()] += conditional_static_cast<N>(p.count);
    }
  });

  return std::make_pair(margs1, margs2);
}
}  // namespace internal

inline std::vector<bool> mad_max_filtering(std::vector<double>& margs, double mad_max) {
  auto mad = [&](const auto vin) {
    const auto median_ = internal::median(vin);
    auto vout = vin;

    std::transform(vout.begin(), vout.end(), vout.begin(),
                   [&](const auto n) { return std::abs(n - median_); });

    return internal::median(vout);
  };

  std::vector<double> cmargs{};
  std::copy_if(margs.begin(), margs.end(), std::back_inserter(cmargs),
               [](const auto n) { return n > 0; });

  if (!cmargs.empty()) {
    const auto median_ = internal::median(cmargs);
    std::transform(margs.begin(), margs.end(), margs.begin(),
                   [&](const auto n) { return n / median_; });
  }

  std::vector<double> log_nz_marg{};
  for (const auto n : margs) {
    if (n > 0) {
      log_nz_marg.push_back(std::log(n));
    }
  }

  std::vector<bool> mask(margs.size(), false);
  if (log_nz_marg.empty()) {
    return mask;
  }

  const auto median_log_nz_marg = internal::median(log_nz_marg);
  const auto dev_log_nz_marg = mad(log_nz_marg);

  const auto cutoff = std::exp(median_log_nz_marg - mad_max * dev_log_nz_marg);

  for (std::size_t i = 0; i < margs.size(); ++i) {
    if (margs[i] < cutoff) {
      mask[i] = true;
    }
  }

  return mask;
}

template <typename PixelIt>
  requires PixelStream<PixelIt>
inline std::pair<std::vector<bool>, std::vector<bool>> mad_max_filtering(
    PixelIt first_pixel, PixelIt last_pixel, const hictk::Chromosome& chrom1,
    const hictk::Chromosome& chrom2, std::uint32_t resolution, double mad_max) {
  if (mad_max == 0) {
    const auto num_bins1 = (chrom1.size() + resolution - 1) / resolution;
    const auto num_bins2 = (chrom2.size() + resolution - 1) / resolution;
    return std::make_pair(std::vector<bool>(num_bins1, false), std::vector<bool>(num_bins2, false));
  }
  auto [margs1, margs2] =
      internal::compute_marginals(first_pixel, last_pixel, chrom1, chrom2, resolution);

  auto mask1 = mad_max_filtering(margs1, mad_max);
  auto mask2 = chrom1 == chrom2 ? mask1 : mad_max_filtering(margs2, mad_max);

  SPDLOG_INFO(
      FMT_STRING("[{}:{}]: MAD-max masking procedure flagged {}/{} bins for {} and {}/{} for {}"),
      chrom1.name(), chrom2.name(), std::accumulate(mask1.begin(), mask1.end(), 0), mask1.size(),
      chrom1.name(), std::accumulate(mask2.begin(), mask2.end(), 0), mask2.size(), chrom2.name());

  return std::make_pair(mask1, mask2);
}

template <typename PixelIt>
  requires PixelStream<PixelIt>
inline std::vector<bool> mad_max_filtering(PixelIt first_pixel, PixelIt last_pixel,
                                           const hictk::Chromosome& chrom, std::uint32_t resolution,
                                           double mad_max) {
  if (mad_max == 0) {
    const auto num_bins = (chrom.size() + resolution - 1) / resolution;
    return std::vector<bool>(num_bins, false);  // NOLINT
  }
  auto [margs, _] = internal::compute_marginals(first_pixel, last_pixel, chrom, chrom, resolution);

  auto mask = mad_max_filtering(margs, mad_max);
  SPDLOG_INFO(FMT_STRING("[{}]: MAD-max masking procedure flagged {}/{} bins"), chrom.name(),
              std::accumulate(mask.begin(), mask.end(), 0), mask.size());
  return mask;
}
}  // namespace nchg
