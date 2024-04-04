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

#include "nchg/expected_matrix.hpp"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>

namespace nchg::test {

TEST_CASE("ExpectedMatrix (cis)", "[short][expected_matrix]") {
  const std::uint32_t resolution = 5;

  const hictk::Chromosome chrom1{0, "chr1", 35};
  const hictk::BinTable bins{{chrom1}, resolution};

  const auto bin1 = bins.at(0);
  const auto bin2 = bins.at(1);
  const auto bin3 = bins.at(2);
  const auto bin4 = bins.at(3);
  const auto bin5 = bins.at(4);
  const auto bin6 = bins.at(5);
  const auto bin7 = bins.at(6);

  using Pixel = hictk::Pixel<std::uint32_t>;
  // clang-format off

  // 1 0 1 0 1 1 1 | 1 0 1 0 1 1 1
  // 0 1 0 0 0 0 0 | 0 1 0 0 0 0 0
  // 0 0 1 0 0 0 0 | 1 0 1 0 0 0 0
  // 0 0 0 1 0 0 0 | 0 0 0 1 0 0 0
  // 0 0 0 0 1 0 0 | 1 0 0 0 1 0 0
  // 0 0 0 0 0 1 0 | 1 0 0 0 0 1 0
  // 0 0 0 0 0 0 1 | 1 0 0 0 0 0 1
  const std::array<Pixel, 11> pixels{
      Pixel{{bin1, bin1}, 1},
      Pixel{{bin1, bin3}, 1},
      Pixel{{bin1, bin5}, 1},
      Pixel{{bin1, bin6}, 1},
      Pixel{{bin1, bin7}, 1},
      Pixel{{bin2, bin2}, 1},
      Pixel{{bin3, bin3}, 1},
      Pixel{{bin4, bin4}, 1},
      Pixel{{bin5, bin5}, 1},
      Pixel{{bin6, bin6}, 1},
      Pixel{{bin7, bin7}, 1}
  };
  // clang-format on

  const ExpectedMatrix m{pixels.begin(), pixels.end(), chrom1, chrom1, bins, 6};

  SECTION("accessors") {
    CHECK(m.resolution() == resolution);
    CHECK(m.num_rows() == 7);
    CHECK(m.num_cols() == 7);
  }

  SECTION("stats") {
    CHECK(m.nnz() == 15);
    CHECK(m.sum() == 15);
    CHECK(m.nnz_avg() == 1);

    REQUIRE(m.marginals1().size() == m.marginals2().size());
    CHECK(std::equal(m.marginals1().begin(), m.marginals1().end(), m.marginals2().begin()));
    CHECK_THAT(m.marginals1().at(0), Catch::Matchers::WithinRel(5.0));
    CHECK_THAT(m.marginals1().at(1), Catch::Matchers::WithinRel(1.0));
    CHECK_THAT(m.marginals1().at(2), Catch::Matchers::WithinRel(2.0));
    CHECK_THAT(m.marginals1().at(3), Catch::Matchers::WithinRel(1.0));
    CHECK_THAT(m.marginals1().at(4), Catch::Matchers::WithinRel(2.0));
    CHECK_THAT(m.marginals1().at(5), Catch::Matchers::WithinRel(2.0));
    CHECK_THAT(m.marginals1().at(6), Catch::Matchers::WithinRel(2.0));
  }

  SECTION("empty matrix") {
    const std::vector<Pixel> buff{};
    const ExpectedMatrix m1{buff.begin(), buff.end(), chrom1, chrom1, bins, 0};
    CHECK(std::all_of(m1.marginals1().begin(), m1.marginals1().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(std::all_of(m1.marginals2().begin(), m1.marginals2().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(m1.nnz() == 0);

    const ExpectedMatrix m2{buff.begin(), buff.end(), chrom1, chrom1, bins, 10};
    CHECK(std::all_of(m2.marginals1().begin(), m2.marginals1().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(std::all_of(m2.marginals2().begin(), m2.marginals2().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(m2.nnz() == 0);
  }
}

TEST_CASE("ExpectedMatrix (trans)", "[short][expected_matrix]") {
  const std::uint32_t resolution = 5;

  const hictk::Chromosome chrom1{0, "chr1", 15};
  const hictk::Chromosome chrom2{1, "chr2", 10};
  const hictk::BinTable bins{{chrom1, chrom2}, resolution};

  const auto bin1 = bins.at(0);
  const auto bin2 = bins.at(1);
  const auto bin3 = bins.at(2);
  const auto bin4 = bins.at(3);
  const auto bin5 = bins.at(4);

  using Pixel = hictk::Pixel<std::uint32_t>;
  // clang-format off

  // 1  0
  // 10 1
  // 0  0
  const std::array<Pixel, 3> pixels{
      Pixel{{bin1, bin4}, 1},
      Pixel{{bin2, bin4}, 10},
      Pixel{{bin2, bin5}, 1},
  };
  // clang-format on

  const ExpectedMatrix m{pixels.begin(), pixels.end(), chrom1, chrom2, bins, 9};

  SECTION("accessors") {
    CHECK(m.resolution() == resolution);
    CHECK(m.num_rows() == 3);
    CHECK(m.num_cols() == 2);

    CHECK(m.at(bin1.id(), bin1.id()) == 4);
    CHECK(m.at(bin1.id(), bin3.id()) == 4);
  }

  SECTION("stats") {
    CHECK(m.nnz() == 3);
    CHECK(m.sum() == 12);
    CHECK(m.nnz_avg() == 4);

    CHECK(m.marginals1().at(0) == 1);
    CHECK(m.marginals1().at(1) == 11);
    CHECK(m.marginals1().at(2) == 0);

    CHECK(m.marginals2().at(0) == 11);
    CHECK(m.marginals2().at(1) == 1);
  }

  SECTION("empty matrix") {
    const std::vector<Pixel> buff{};
    const ExpectedMatrix m1{buff.begin(), buff.end(), chrom1, chrom1, bins, 0};
    CHECK(std::all_of(m1.marginals1().begin(), m1.marginals1().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(std::all_of(m1.marginals2().begin(), m1.marginals2().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(m1.nnz() == 0);

    const ExpectedMatrix m2{buff.begin(), buff.end(), chrom1, chrom1, bins, 10};
    CHECK(std::all_of(m2.marginals1().begin(), m2.marginals1().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(std::all_of(m2.marginals2().begin(), m2.marginals2().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(m2.nnz() == 0);
  }
}

}  // namespace nchg::test
