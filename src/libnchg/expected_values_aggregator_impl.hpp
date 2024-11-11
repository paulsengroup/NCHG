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

#include <cmath>
#include <hictk/pixel.hpp>

#include "nchg/common.hpp"

namespace nchg {

template <typename N>
inline void ExpectedValuesAggregator::add(const hictk::ThinPixel<N> &p) {
  add(hictk::Pixel<N>{*_bins, p});
}

template <typename N>
inline void ExpectedValuesAggregator::add(const hictk::Pixel<N> &p) {
  const auto count = conditional_static_cast<double>(p.count);
  if (std::isnan(count)) [[unlikely]] {
    return;
  }

  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();

  if (p.coords.is_intra()) {
    at(chrom1) += count;
    const auto i = p.coords.bin2.id() - p.coords.bin1.id();
    // skip last bin in chromosome if chromosome size is not a multiple of bin size
    // this is done to mimic HiCTools' behavior
    if (i < _observed_distances.size()) [[likely]] {
      _observed_distances[i] += count;
    }
  } else {
    at(chrom1, chrom2) += count;
  }
}

}  // namespace nchg
