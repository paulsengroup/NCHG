
// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
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

#include "nchg/genomic_domains.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/pixel.hpp>

namespace nchg::test {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("GenomicDomains", "[short][genomic_domains]") {
  const hictk::Chromosome chrom1{0, "chr1", 1'000};
  const hictk::Chromosome chrom2{1, "chr2", 500};
  const hictk::Chromosome chrom3{2, "chr3", 100};

  const hictk::GenomicInterval gi1{chrom1, 10, 20};
  const hictk::GenomicInterval gi2{chrom2, 5, 15};

  const GenomicDomains domains{{
      {gi1, gi1},
      {gi1, gi2},
      {gi2, gi2},
  }};

  SECTION("ctor w/ duplicates") {
    CHECK(GenomicDomains{{
                             {gi1, gi1},
                             {gi1, gi1},
                             {gi2, gi2},
                         }}
              .size() == 2);
  }

  SECTION("accessors") {
    CHECK(GenomicDomains{}.empty());
    CHECK(domains.size() == 3);
    CHECK(domains().size() == 3);

    CHECK(domains.contains(chrom1));
    CHECK(domains.contains(chrom2));
    CHECK(domains.contains(chrom1, chrom2));
    CHECK(!domains.contains(chrom3));
    CHECK(!domains.contains(chrom1, chrom3));
  }

  SECTION("fetch") {
    CHECK(domains.fetch<int>(chrom1).size() == 1);
    CHECK(domains.fetch<int>(chrom1, chrom2).size() == 1);
    CHECK(domains.fetch<int>(chrom3).empty());
  }
}

TEST_CASE("GenomicDomainsIndexed", "[short][genomic_domains]") {
  const hictk::Chromosome chrom1{0, "chr1", 1'000};

  SECTION("1D") {
    auto make_domain = [&](std::uint32_t start, std::uint32_t end) {
      return BEDPE{{chrom1, start, end}, {chrom1, start, end}};
    };

    auto make_pixel = [&](std::uint32_t start, std::uint32_t end, int count = 1) {
      return hictk::Pixel{chrom1, start, end, count};
    };

    const hictk::GenomicInterval gi1{chrom1, 10, 20};
    const hictk::GenomicInterval gi2{chrom1, 15, 25};
    const hictk::GenomicInterval gi3{chrom1, 20, 30};
    const hictk::GenomicInterval gi4{chrom1, 50, 60};

    auto domains = GenomicDomains{
        {
            {gi1, gi1},
            {gi2, gi2},
            {gi3, gi3},
            {gi4, gi4},
        }}.fetch<int>(chrom1);

    REQUIRE(domains.size() == 4);

    SECTION("empty") {
      for (const auto& [_, n] : domains.to_vector()) {
        CHECK(n == 0);
      }
    }

    SECTION("overlaps") {
      domains.add_interactions(make_pixel(50, 60));
      CHECK(domains.sum() == 1);
      CHECK(domains.at(make_domain(50, 60)) == 1);
    }

    SECTION("no overlap") {
      domains.add_interactions(make_pixel(60, 70));
      CHECK(domains.sum() == 0);
    }

    SECTION("partial overlap") {
      domains.add_interactions(make_pixel(10, 12));
      CHECK(domains.sum() == 1);
      CHECK(domains.at(make_domain(10, 20)) == 1);
    }

    SECTION("multiple overlaps") {
      domains.add_interactions(make_pixel(15, 20));
      CHECK(domains.sum() == 2);
      CHECK(domains.at(make_domain(10, 20)) == 1);
      CHECK(domains.at(make_domain(15, 25)) == 1);
    }
  }

  SECTION("2D") {
    auto make_domain = [&](std::uint32_t start1, std::uint32_t end1, std::uint32_t start2,
                           std::uint32_t end2) {
      return BEDPE{{chrom1, start1, end1}, {chrom1, start2, end2}};
    };

    auto make_pixel = [&](std::uint32_t start1, std::uint32_t end1, std::uint32_t start2,
                          std::uint32_t end2, int count = 1) {
      return hictk::Pixel{chrom1, start1, end1, chrom1, start2, end2, count};
    };

    const hictk::GenomicInterval gi1{chrom1, 0, 10};
    const hictk::GenomicInterval gi2{chrom1, 10, 20};
    const hictk::GenomicInterval gi3{chrom1, 15, 25};
    const hictk::GenomicInterval gi4{chrom1, 16, 18};
    const hictk::GenomicInterval gi5{chrom1, 50, 60};

    auto domains = GenomicDomains{
        {
            {gi1, gi5},
            {gi2, gi5},
            {gi3, gi5},
            {gi4, gi5},
        }}.fetch<int>(chrom1);

    REQUIRE(domains.size() == 4);

    SECTION("empty") {
      for (const auto& [_, n] : domains.to_vector()) {
        CHECK(n == 0);
      }
    }

    SECTION("overlaps") {
      domains.add_interactions(make_pixel(0, 10, 50, 60));
      CHECK(domains.sum() == 1);
      CHECK(domains.at(make_domain(0, 10, 50, 60)) == 1);
    }

    SECTION("no overlap") {
      domains.add_interactions(make_pixel(10, 20, 10, 20));
      CHECK(domains.sum() == 0);

      domains.add_interactions(make_pixel(50, 60, 50, 60));
      CHECK(domains.sum() == 0);
    }

    SECTION("partial overlap") {
      SECTION("first dimension") {
        domains.add_interactions(make_pixel(10, 15, 50, 60));
        CHECK(domains.sum() == 1);
        CHECK(domains.at(make_domain(10, 20, 50, 60)) == 1);
      }

      SECTION("second dimension") {
        domains.add_interactions(make_pixel(0, 10, 55, 60));
        CHECK(domains.sum() == 1);
        CHECK(domains.at(make_domain(0, 10, 50, 60)) == 1);
      }

      SECTION("both dimensions") {
        domains.add_interactions(make_pixel(0, 5, 55, 60));
        CHECK(domains.sum() == 1);
        CHECK(domains.at(make_domain(0, 10, 50, 60)) == 1);
      }
    }

    SECTION("multiple overlaps") {
      domains.add_interactions(make_pixel(10, 20, 50, 60));
      CHECK(domains.sum() == 3);
      CHECK(domains.at(make_domain(10, 20, 50, 60)) == 1);
      CHECK(domains.at(make_domain(15, 25, 50, 60)) == 1);
      CHECK(domains.at(make_domain(16, 18, 50, 60)) == 1);
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace nchg::test
