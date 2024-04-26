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

#include <spdlog/spdlog.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <hictk/cooler/cooler.hpp>
#include <hictk/fmt/pixel.hpp>
#include <memory>

#include "tmpdir.hpp"

namespace nchg::test {

TEST_CASE("NCHG", "[long][nchg]") {
  const auto test_file = datadir / "ENCFF447ERX.1000000.cool";

  const auto clr = std::make_shared<const hictk::cooler::File>(test_file.string());

  auto params = NCHG<hictk::cooler::File>::DefaultParams;
  params.mad_max = 0.0;

  const auto chr1 = clr->chromosomes().at("chr1");
  const NCHG nchg(clr, chr1, chr1, params);

  SECTION("significant") {
    const hictk::GenomicInterval range1(chr1, 0, 10'000'000);
    const hictk::GenomicInterval range2(chr1, 10'000'000, 20'000'000);

    const auto s = nchg.compute(range1, range2, 1.0);

    CHECK(s.pval <= 0.05);
  }
  SECTION("not significant") {
    const hictk::GenomicInterval range1(chr1, 0, 10'000'000);
    const hictk::GenomicInterval range2(chr1, 95'000'000, 105'000'000);

    const auto s = nchg.compute(range1, range2, 1.0);

    CHECK(s.pval > 0.05);
  }
}

}  // namespace nchg::test
