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

#include "nchg/expected_values.hpp"

#include <spdlog/spdlog.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <hictk/file.hpp>
#include <memory>

#include "nchg/test/tmpdir.hpp"

namespace nchg::test {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
static void compare_weights(const std::vector<double>& weights, const std::vector<double>& expected,
                            double tol = 1.0e-6) {
  REQUIRE(weights.size() == expected.size());

  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (std::isnan(weights[i])) {
      CHECK(std::isnan(expected[i]));
    } else {
      CHECK_THAT(weights[i], Catch::Matchers::WithinRel(expected[i], tol));
    }
  }
}

static void compare_masks(const std::vector<bool>& mask, const std::vector<bool>& expected) {
  REQUIRE(mask.size() == expected.size());

  for (std::size_t i = 0; i < mask.size(); ++i) {
    CHECK(mask[i] == expected[i]);
  }
}

[[nodiscard]] static std::vector<double> get_weights(
    const ExpectedValues& evs, const std::shared_ptr<const hictk::File>& clr) {
  return evs.expected_values(clr->chromosomes().longest_chromosome(), false);
}

TEST_CASE("ExpectedValues: genome-wide", "[long][expected_values]") {
  spdlog::default_logger()->set_level(spdlog::level::warn);
  const auto test_file = datadir / "ENCFF447ERX.1000000.cool";

  const auto clr = std::make_shared<const hictk::File>(test_file.string());

  const ExpectedValues evs(clr);

  SECTION("serde") {
    const auto path = testdir() / "expected_values_gw.bin";
    evs.serialize(path);
    const auto evs_gw_serde = ExpectedValues::deserialize(path);

    CHECK(evs.params().mad_max == evs_gw_serde.params().mad_max);
    CHECK(evs.params().min_delta == evs_gw_serde.params().min_delta);
    CHECK(evs.params().max_delta == evs_gw_serde.params().max_delta);
    CHECK(evs.params().bin_aggregation_observed_distances_cutoff ==
          evs_gw_serde.params().bin_aggregation_observed_distances_cutoff);
    CHECK(evs.params().bin_aggregation_possible_distances_cutoff ==
          evs_gw_serde.params().bin_aggregation_possible_distances_cutoff);
    CHECK(evs.params().interpolate == evs_gw_serde.params().interpolate);
    CHECK(evs.params().interpolation_qtile == evs_gw_serde.params().interpolation_qtile);
    CHECK(evs.params().interpolation_window_size ==
          evs_gw_serde.params().interpolation_window_size);

    const auto w1 = get_weights(evs, clr);
    const auto w2 = get_weights(evs_gw_serde, clr);

    compare_weights(w1, w2);

    CHECK(evs.scaling_factors().size() == evs_gw_serde.scaling_factors().size());
    for (const auto& [chrom, sf] : evs.scaling_factors()) {
      CHECK(evs_gw_serde.scaling_factors().at(chrom) == sf);
    }

    const auto chrom1 = clr->chromosomes().at("chr21");
    const auto chrom2 = clr->chromosomes().at("chr22");
    CHECK(evs.expected_value(chrom1, chrom2) == evs_gw_serde.expected_value(chrom1, chrom2));

    compare_masks(*evs_gw_serde.bin_mask(chrom1, chrom2).first,
                  *evs.bin_mask(chrom1, chrom2).first);
    compare_masks(*evs_gw_serde.bin_mask(chrom1, chrom2).second,
                  *evs.bin_mask(chrom1, chrom2).second);
  }
}

TEST_CASE("ExpectedValues: cis-only", "[medium][expected_values]") {
  spdlog::default_logger()->set_level(spdlog::level::warn);
  const auto test_file = datadir / "ENCFF447ERX.1000000.cool";

  const auto clr = std::make_shared<const hictk::File>(test_file.string());

  const auto evs = ExpectedValues::cis_only(clr);

  SECTION("accessors") {
    const auto num_bins =
        (clr->chromosomes().at("chr1").size() + clr->resolution() - 1) / clr->resolution();
    CHECK(get_weights(evs, clr).size() == num_bins);

    const auto chrom = clr->chromosomes().at("chr22");
    const auto chrom_num_bins = (chrom.size() + clr->resolution() - 1) / clr->resolution();
    CHECK(evs.expected_values(chrom).size() == chrom_num_bins);

    CHECK(evs.scaling_factors().size() == clr->chromosomes().size());

    CHECK_THROWS(evs.expected_value(chrom, chrom));
  }

  SECTION("expected values") {
    const auto chrom1 = clr->chromosomes().at("chr21");
    const auto chrom2 = clr->chromosomes().at("chr22");

    CHECK_NOTHROW(evs.expected_values(chrom1));
    CHECK_THROWS(evs.expected_value(chrom1, chrom2));

    const auto num_bins = (chrom1.size() + clr->resolution() - 1) / clr->resolution();
    CHECK(evs.expected_matrix(chrom1).num_rows() == num_bins);
  }
}

TEST_CASE("ExpectedValues: trans-only", "[long][expected_values]") {
  spdlog::default_logger()->set_level(spdlog::level::warn);
  const auto test_file = datadir / "ENCFF447ERX.1000000.cool";

  const auto clr = std::make_shared<const hictk::File>(test_file.string());

  auto params = ExpectedValues::DefaultParams;
  params.mad_max = 0.0;

  const auto evs = ExpectedValues::trans_only(clr, params);
  const auto chrom1 = clr->chromosomes().at("chr21");
  const auto chrom2 = clr->chromosomes().at("chr22");

  CHECK_THROWS(evs.expected_values(chrom1));
  CHECK_NOTHROW(evs.expected_value(chrom1, chrom2));
}

TEST_CASE("ExpectedValues: chromosome pair", "[medium][expected_values]") {
  spdlog::default_logger()->set_level(spdlog::level::warn);
  const auto test_file = datadir / "ENCFF447ERX.1000000.cool";

  const auto clr = std::make_shared<const hictk::File>(test_file.string());
  const auto chrom1 = clr->chromosomes().at("chr21");
  const auto chrom2 = clr->chromosomes().at("chr22");
  const auto chrom3 = clr->chromosomes().at("chrX");

  const auto evs1 = ExpectedValues::chromosome_pair(clr, chrom1, chrom1);
  const auto evs2 = ExpectedValues::chromosome_pair(clr, chrom1, chrom2);

  SECTION("expected values") {
    CHECK_NOTHROW(evs1.expected_values(chrom1));
    CHECK_THROWS(evs2.expected_values(chrom1));

    CHECK_THROWS(evs1.expected_value(chrom1, chrom2));
    CHECK_NOTHROW(evs2.expected_value(chrom1, chrom2));
    CHECK_THROWS(evs2.expected_value(chrom1, chrom3));
  }

  SECTION("expected matrix") {
    CHECK(evs2.expected_matrix(chrom1, chrom2).nnz_avg() == evs2.expected_value(chrom1, chrom2));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace nchg::test
