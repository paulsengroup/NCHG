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

#include "nchg/observed_matrix.hpp"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/pixel.hpp>

namespace nchg::test {

TEST_CASE("ObservedMatrix (cis)", "[short][observed_matrix]") {
  const std::uint32_t resolution = 5;

  const hictk::Chromosome chrom1{0, "chr1", 15};
  const hictk::BinTable bins{{chrom1}, resolution};

  const auto bin1 = bins.at(0);
  const auto bin2 = bins.at(1);
  const auto bin3 = bins.at(2);

  using Pixel = hictk::Pixel<std::uint32_t>;
  // clang-format off

  // 1 0 1 | 1 0 1
  // 0 1 0 | 0 1 0
  // 0 0 1 | 1 0 1
  const std::array<Pixel, 4> pixels{
      Pixel{{bin1, bin1}, 1},
      Pixel{{bin1, bin3}, 1},
      Pixel{{bin2, bin2}, 1},
      Pixel{{bin3, bin3}, 1}
  };
  // clang-format on

  SECTION("full matrix") {
    const ObservedMatrix m{pixels.begin(), pixels.end(), chrom1, chrom1, bins};

    SECTION("accessors") {
      CHECK(m.resolution() == resolution);
      CHECK(m.num_rows() == 3);
      CHECK(m.num_cols() == 3);
    }

    SECTION("stats") {
      CHECK(m.nnz() == 5);
      CHECK(m.sum() == 5);
      CHECK(m.nnz_avg() == 1);

      REQUIRE(m.marginals1().size() == m.marginals2().size());
      CHECK(std::equal(m.marginals1().begin(), m.marginals1().end(), m.marginals2().begin()));
      CHECK(m.marginals1().at(0) == 2);
      CHECK(m.marginals1().at(1) == 1);
      CHECK(m.marginals1().at(2) == 2);
    }

    SECTION("iterators") {
      const auto nnz = static_cast<std::ptrdiff_t>(pixels.size());
      CHECK(std::distance(m.begin(), m.end()) == nnz);
    }
  }

  SECTION("masked matrix") {
    const ObservedMatrix m{pixels.begin(), pixels.end(), chrom1, chrom1, bins, 1, 100};

    SECTION("accessors") {
      CHECK(m.resolution() == resolution);
      CHECK(m.num_rows() == 3);
      CHECK(m.num_cols() == 3);

      CHECK(m.min_delta() == 1);
      CHECK(m.max_delta() == 100);
    }

    SECTION("stats") {
      CHECK(m.nnz() == 2);
      CHECK(m.sum() == 2);
      CHECK(m.nnz_avg() == 1);

      REQUIRE(m.marginals1().size() == m.marginals2().size());
      CHECK(std::equal(m.marginals1().begin(), m.marginals1().end(), m.marginals2().begin()));
      CHECK(m.marginals1().at(0) == 1);
      CHECK(m.marginals1().at(1) == 0);
      CHECK(m.marginals1().at(2) == 1);
    }

    SECTION("iterators") {
      const auto nnz = static_cast<std::ptrdiff_t>(pixels.size());
      CHECK(std::distance(m.begin(), m.end()) == nnz);
    }
  }

  SECTION("empty matrix") {
    const std::vector<Pixel> buff{};
    const ObservedMatrix m{buff.begin(), buff.end(), chrom1, chrom1, bins};
    CHECK(std::all_of(m.marginals1().begin(), m.marginals1().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(std::all_of(m.marginals2().begin(), m.marginals2().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(m.nnz() == 0);
  }
}

TEST_CASE("ObservedMatrix (trans)", "[short][observed_matrix]") {
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

  // 1 0
  // 1 1
  // 0 0
  const std::array<Pixel, 3> pixels{
      Pixel{{bin1, bin4}, 1},
      Pixel{{bin2, bin4}, 1},
      Pixel{{bin2, bin5}, 1},
  };
  // clang-format on

  SECTION("full matrix") {
    const ObservedMatrix m{pixels.begin(), pixels.end(), chrom1, chrom2, bins};

    SECTION("accessors") {
      CHECK(m.resolution() == resolution);
      CHECK(m.num_rows() == 3);
      CHECK(m.num_cols() == 2);
    }

    SECTION("stats") {
      CHECK(m.nnz() == 3);
      CHECK(m.sum() == 3);
      CHECK(m.nnz_avg() == 1);

      CHECK(m.marginals1().at(0) == 1);
      CHECK(m.marginals1().at(1) == 2);
      CHECK(m.marginals1().at(2) == 0);

      CHECK(m.marginals2().at(0) == 2);
      CHECK(m.marginals2().at(1) == 1);
    }

    SECTION("iterators") {
      const auto nnz = static_cast<std::ptrdiff_t>(pixels.size());
      CHECK(std::distance(m.begin(), m.end()) == nnz);
    }
  }

  SECTION("masked matrix") {
    const ObservedMatrix m{pixels.begin(), pixels.end(), chrom1, chrom2, bins, 1, 100};

    SECTION("accessors") {
      CHECK(m.resolution() == resolution);
      CHECK(m.num_rows() == 3);
      CHECK(m.num_cols() == 2);

      CHECK(m.min_delta() == 1);
      CHECK(m.max_delta() == 100);
    }

    SECTION("stats") {
      CHECK(m.nnz() == 3);
      CHECK(m.sum() == 3);
      CHECK(m.nnz_avg() == 1);

      CHECK(m.marginals1().at(0) == 1);
      CHECK(m.marginals1().at(1) == 2);
      CHECK(m.marginals1().at(2) == 0);

      CHECK(m.marginals2().at(0) == 2);
      CHECK(m.marginals2().at(1) == 1);
    }

    SECTION("iterators") {
      const auto nnz = static_cast<std::ptrdiff_t>(pixels.size());
      CHECK(std::distance(m.begin(), m.end()) == nnz);
    }

    SECTION("empty matrix") {
      const std::vector<Pixel> buff{};
      const ObservedMatrix m1{buff.begin(), buff.end(), chrom1, chrom1, bins};
      CHECK(std::all_of(m1.marginals1().begin(), m1.marginals1().end(),
                        [](const auto n) { return n == 0; }));
      CHECK(std::all_of(m1.marginals2().begin(), m1.marginals2().end(),
                        [](const auto n) { return n == 0; }));
      CHECK(m1.nnz() == 0);
    }
  }

  SECTION("empty matrix") {
    const std::vector<Pixel> buff{};
    const ObservedMatrix m{buff.begin(), buff.end(), chrom1, chrom1, bins};
    CHECK(std::all_of(m.marginals1().begin(), m.marginals1().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(std::all_of(m.marginals2().begin(), m.marginals2().end(),
                      [](const auto n) { return n == 0; }));
    CHECK(m.nnz() == 0);
  }
}

}  // namespace nchg::test
