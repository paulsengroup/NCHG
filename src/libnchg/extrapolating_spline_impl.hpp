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

#include <algorithm>
#include <boost/math/interpolators/pchip.hpp>
#include <cassert>
#include <cstdint>
#include <vector>

// NOLINTNEXTLINE(*-concat-nested-namespaces)
namespace nchg {

namespace internal {
[[nodiscard]] inline std::pair<double, double> extrapolate_x_intercept(const std::vector<double>& x,
                                                                       const std::vector<double>& y,
                                                                       std::size_t n = 1) {
  assert(x.size() == y.size());
  assert(x.size() > n);

  const auto x1 = x.front();
  const auto x2 = x.at(n);
  const auto y1 = y.front();
  const auto y2 = y.at(n);

  const auto m = (y1 - y2) / (x1 - x2);

  const auto x0 = 0.0;
  const auto y0 = y1 + (m * (x0 - x1));
  return std::make_pair(x0, y0);
}

[[nodiscard]] inline std::pair<double, double> extrapolate_y_intercept(const std::vector<double>& x,
                                                                       const std::vector<double>& y,
                                                                       std::size_t n = 1) {
  assert(x.size() == y.size());
  assert(x.size() > n);

  const auto x1 = x.at(x.size() - (n + 1));
  const auto x2 = x.back();
  const auto y1 = y.at(y.size() - (n + 1));
  const auto y2 = y.back();

  const auto m = (y2 - y1) / (x2 - x1);

  const auto yn = 0.0;
  const auto xn = x2 + ((yn - y2) / m);
  return std::make_pair(xn, yn);
}

[[nodiscard]] inline std::pair<double, double> extrapolate_y(const std::vector<double>& x,
                                                             const std::vector<double>& y,
                                                             double xn, std::size_t n = 1) {
  assert(x.size() == y.size());
  assert(x.size() > n);
  assert(xn >= x.back());

  const auto x1 = x.front();
  const auto x2 = x.at(n);
  const auto y1 = y.front();
  const auto y2 = y.at(n);

  const auto m = (y1 - y2) / (x1 - x2);

  const auto y0 = y1 + (m * (xn - x1));
  return std::make_pair(xn, y0);
}

}  // namespace internal

inline Extrapolating1DSpline::Extrapolating1DSpline(std::vector<double> x, std::vector<double> y) {
  if (x.front() != 0) {
    const auto [x0, y0] = internal::extrapolate_x_intercept(x, y);
    if (std::isfinite(x0)) {
      x.insert(x.begin(), x0);
      y.insert(y.begin(), y0);
    }
  }

  if (y.back() > 0) {
    double xn = -1;
    double yn = -1;

    for (std::size_t i = 1; i < x.size() - 1 && xn < x.back(); ++i) {
      const auto res = internal::extrapolate_y_intercept(x, y, i);
      xn = res.first;
      yn = res.second;
    }
    if (xn >= x.back()) {
      x.push_back(xn);
      y.push_back(yn);
    } else {
      // y should be monotonically decreasing, but just in case it isn't we extrapolate y for 2 *
      // x.back()
      const auto [xn_, yn_] = internal::extrapolate_y(x, y, 2 * x.back());
      x.push_back(xn_);
      y.push_back(yn_);
    }
  }

  _pchip = boost::math::interpolators::pchip(std::move(x), std::move(y));
}

[[nodiscard]] inline double Extrapolating1DSpline::evaluate(std::uint64_t x) const {
  const auto xd = static_cast<double>(x);
  return _pchip.value()(xd);
}

}  // namespace nchg
