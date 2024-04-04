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

#include "nchg/extrapolating_spline.hpp"

#include <catch2/catch_test_macros.hpp>
#include <vector>

namespace nchg::test {

TEST_CASE("ExtrapolatingSpline", "[short][extrapolating_spline]") {
  const std::vector<double> x{2, 4, 6, 8, 10};
  const std::vector<double> y = {10, 8, 6, 4, 2};  // NOLINT

  const Extrapolating1DSpline spline{x, y};

  CHECK(spline.evaluate(0) == 12);
  CHECK(spline.evaluate(1) == 11);
  CHECK(spline.evaluate(2) == 10);
  CHECK(spline.evaluate(11) == 1);
}

}  // namespace nchg::test
