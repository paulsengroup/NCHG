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

#include <cstdint>
#include <hictk/chromosome.hpp>
#include <vector>

#include "nchg/concepts.hpp"

namespace nchg {

[[nodiscard]] std::vector<bool> mad_max_filtering(std::vector<double>& margs, double mad_max);

template <typename PixelIt>
  requires PixelStream<PixelIt>
[[nodiscard]] std::vector<bool> mad_max_filtering(PixelIt first_pixel, PixelIt last_pixel,
                                                  const hictk::Chromosome& chrom,
                                                  std::uint32_t resolution, double mad_max);

template <typename PixelIt>
  requires PixelStream<PixelIt>
[[nodiscard]] std::pair<std::vector<bool>, std::vector<bool>> mad_max_filtering(
    PixelIt first_pixel, PixelIt last_pixel, const hictk::Chromosome& chrom1,
    const hictk::Chromosome& chrom2, std::uint32_t resolution, double mad_max);

}  // namespace nchg

#include "../../mad_max_filter_impl.hpp"
