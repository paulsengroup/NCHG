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

#include "nchg/mad_max_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <vector>

#include "nchg/median.hpp"

namespace nchg {

std::vector<bool> mad_max_filtering(std::vector<double>& margs, double mad_max) {
  auto mad = [&](const auto& vin) {
    const auto median_ = median(vin);
    auto vout = vin;

    std::ranges::transform(vout, vout.begin(), [&](const auto n) { return std::abs(n - median_); });

    return median(vout);
  };

  if (std::ranges::all_of(margs, [](const auto n) { return n == 0; })) {
    return std::vector(margs.size(), true);
  }

  std::vector<double> cmargs{};
  std::ranges::copy_if(margs, std::back_inserter(cmargs), [](const auto n) { return n > 0; });

  if (!cmargs.empty()) {
    const auto median_ = median(cmargs);
    std::ranges::transform(margs, margs.begin(), [&](const auto n) { return n / median_; });
  }

  std::vector<double> log_nz_marg{};
  for (const auto n : margs) {
    if (n > 0) {
      log_nz_marg.push_back(std::log(n));
    }
  }

  std::vector mask(margs.size(), false);
  if (log_nz_marg.empty()) {
    return mask;
  }

  const auto median_log_nz_marg = median(log_nz_marg);
  const auto dev_log_nz_marg = mad(log_nz_marg);

  const auto cutoff = std::exp(median_log_nz_marg - (mad_max * dev_log_nz_marg));

  for (std::size_t i = 0; i < margs.size(); ++i) {
    if (margs[i] < cutoff) {
      mask[i] = true;
    }
  }

  return mask;
}

}  // namespace nchg
