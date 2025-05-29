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

#include "nchg/nchg.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <hictk/file.hpp>
#include <memory>

#include "nchg/genomic_domains.hpp"
#include "nchg/test/tmpdir.hpp"

namespace nchg::test {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("NCHG", "[medium][nchg]") {
  const auto test_file = datadir / "ENCFF447ERX.minified.mcool";

  const auto clr = std::make_shared<const hictk::File>(test_file.string(), 1'000'000);

  auto params = NCHG::DefaultParams;
  params.mad_max = 0.0;

  SECTION("cis") {
    const auto chr21 = clr->chromosomes().at("chr21");
    const NCHG nchg(clr, chr21, chr21, params);

    SECTION("significant") {
      const BEDPE domain{{chr21, 6'000'000, 7'000'000}, {chr21, 43'000'000, 44'000'000}};

      constexpr std::uint64_t obs = 335769;
      constexpr double exp = 12578.709427;
      const auto s = nchg.compute(domain, obs, exp, 1.0);

      CHECK(s.pval <= 0.05);
      CHECK(s.log_ratio > 2.0);
    }
    SECTION("not significant") {
      const BEDPE domain{{chr21, 6'000'000, 7'000'000}, {chr21, 10'000'000, 11'000'000}};

      constexpr std::uint64_t obs = 1945;
      constexpr double exp = 19098.256428;
      const auto s = nchg.compute(domain, obs, exp, 1.0);

      CHECK(s.pval > 0.05);
      CHECK(s.log_ratio < 2.0);
    }
  }

  SECTION("trans") {
    const auto chr21 = clr->chromosomes().at("chr21");
    const auto chr22 = clr->chromosomes().at("chr22");
    const NCHG nchg(clr, chr21, chr22, params);

    SECTION("significant") {
      const BEDPE domain{{chr21, 10'000'000, 11'000'000}, {chr22, 11'000'000, 12'000'000}};

      constexpr std::uint64_t obs = 165978;
      constexpr double exp = 1859.162821;
      const auto s = nchg.compute(domain, obs, exp, 1.0);

      CHECK(s.pval <= 0.05);
      CHECK(s.log_ratio > 2.0);
      CHECK(s.omega == 1);
    }
    SECTION("not significant") {
      const BEDPE domain{{chr21, 21'000'000, 22'000'000}, {chr22, 22'000'000, 23'000'000}};

      constexpr std::uint64_t obs = 163;
      constexpr double exp = 1859.162821;
      const auto s = nchg.compute(domain, obs, exp, 1.0);

      CHECK(s.pval > 0.05);
      CHECK(s.log_ratio < 2.0);
      CHECK(s.omega == 1);
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace nchg::test
