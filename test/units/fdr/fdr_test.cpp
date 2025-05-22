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

#include "nchg/fdr.hpp"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>

namespace nchg::test {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("FDR", "[short][fdr]") {
  SECTION("floating point") {
    SECTION("unique values") {
      const std::vector<double> pvalues{0.2500000000000000, 0.0000000000000000, 0.0833333333333333,
                                        0.1666666666666667, 0.0555555555555556, 0.1111111111111111,
                                        0.1388888888888889, 0.0277777777777778, 0.1944444444444444,
                                        0.2222222222222222};
      const std::vector<double> corrected_pvalues_expected{
          0.2500000000000000, 0.0000000000000000, 0.2083333333333333, 0.2380952380952381,
          0.1851851851851852, 0.2222222222222222, 0.2314814814814815, 0.1388888888888889,
          0.2430555555555555, 0.2469135802469136};
      BH_FDR fdr(pvalues);

      const auto corrected_pvalues = fdr.correct();
      REQUIRE(corrected_pvalues_expected.size() == corrected_pvalues.size());
      for (std::size_t i = 0; i < corrected_pvalues.size(); ++i) {
        CHECK_THAT(corrected_pvalues[i], Catch::Matchers::WithinRel(corrected_pvalues_expected[i]));
      }
    }

    SECTION("duplicate values") {
      const std::vector<double> pvalues{0.05, 0.05, 0.1, 0.0};
      const std::vector<double> corrected_pvalues_expected{0.0666666666666667, 0.0666666666666667,
                                                           0.1, 0.0};
      BH_FDR fdr(pvalues);

      const auto corrected_pvalues = fdr.correct();
      REQUIRE(corrected_pvalues_expected.size() == corrected_pvalues.size());
      for (std::size_t i = 0; i < corrected_pvalues.size(); ++i) {
        CHECK_THAT(corrected_pvalues[i], Catch::Matchers::WithinRel(corrected_pvalues_expected[i]));
      }
    }
  }

  SECTION("struct") {
    struct Stats {
      std::size_t foo{};
      double pval{};
    };

    // NOLINTNEXTLINE(*-use-designated-initializers)
    const std::vector<Stats> pvalues{{0, 0.05}, {1, 0.1}, {2, 0.0}};
    const std::vector corrected_pvalues_expected{0.075, 0.1, 0.0};

    BH_FDR fdr(pvalues);

    const auto corrected_pvalues = fdr.correct([](Stats& s) -> double& { return s.pval; });
    REQUIRE(corrected_pvalues_expected.size() == corrected_pvalues.size());
    for (std::size_t i = 0; i < corrected_pvalues.size(); ++i) {
      CHECK(corrected_pvalues[i].foo == i);
      CHECK_THAT(corrected_pvalues[i].pval,
                 Catch::Matchers::WithinRel(corrected_pvalues_expected[i]));
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace nchg::test
