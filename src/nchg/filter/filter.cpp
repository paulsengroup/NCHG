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

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "nchg/fdr.hpp"
#include "nchg/tools.hpp"

namespace nchg {

static void print_header() {
  fmt::print(
      FMT_STRING("chrom1\t"
                 "start1\t"
                 "end1\t"
                 "chrom2\t"
                 "start2\t"
                 "end2\t"
                 "pvalue\t"
                 "pvalue_corrected\t"
                 "observed_count\t"
                 "expected_count\t"
                 "log_ratio\t"
                 "odds_ratio\t"
                 "omega\n"));
}

struct SharedPtrStringCmp {
  using is_transparent = void;
  bool operator()(const std::shared_ptr<std::string>& a,
                  const std::shared_ptr<std::string>& b) const noexcept {
    assert(!!a);
    assert(!!b);
    return *a < *b;
  }

  bool operator()(const std::string& a, const std::shared_ptr<std::string>& b) const noexcept {
    assert(!!b);
    return a < *b;
  }

  bool operator()(const std::shared_ptr<std::string>& a, const std::string& b) const noexcept {
    assert(!!a);
    return *a < b;
  }

  bool operator()(std::string_view a, const std::shared_ptr<std::string>& b) const noexcept {
    assert(!!b);
    return a < *b;
  }

  bool operator()(const std::shared_ptr<std::string>& a, std::string_view b) const noexcept {
    assert(!!a);
    return *a < b;
  }
};

using ChromosomeSet = phmap::btree_set<std::shared_ptr<std::string>, SharedPtrStringCmp>;

struct Stats {
  std::shared_ptr<const std::string> chrom1{};
  std::uint32_t start1{};
  std::uint32_t end1{};

  std::shared_ptr<const std::string> chrom2{};
  std::uint32_t start2{};
  std::uint32_t end2{};

  double pvalue{};
  double pvalue_corrected{};
  std::uint64_t observed_count{};
  double expected_count{};
  double log_ratio{};
  double odds_ratio{};
  double omega{};

  [[nodiscard]] static Stats parse(ChromosomeSet& chromosomes, std::string_view s,
                                   char sep = '\t') {
    assert(!s.empty());
    assert(s.find('#') != 0);
    // TODO make more efficient
    std::vector<std::string_view> toks{};
    while (!s.empty()) {
      const auto pos = s.find(sep);
      toks.push_back(s.substr(0, pos));
      if (pos == std::string_view::npos) {
        break;
      }
      s = s.substr(pos + 1);
    }

    if (toks.size() != 11) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected 11 fields, found {}"), toks.size()));
    }

    if (chromosomes.find(toks[0]) == chromosomes.end()) {
      chromosomes.emplace(std::make_shared<std::string>(std::string{toks[0]}));
    }
    auto chrom1_ptr = *chromosomes.find(toks[0]);
    const auto start1 = std::stoull(std::string{toks[1]});
    const auto end1 = std::stoull(std::string{toks[2]});

    if (chromosomes.find(toks[3]) == chromosomes.end()) {
      chromosomes.emplace(std::make_shared<std::string>(std::string{toks[3]}));
    }
    auto chrom2_ptr = *chromosomes.find(toks[3]);
    const auto start2 = std::stoull(std::string{toks[4]});
    const auto end2 = std::stoull(std::string{toks[5]});

    const auto pvalue = std::stod(std::string{toks[6]});
    const auto obs = std::stoull(std::string{toks[7]});
    const auto exp = std::stod(std::string{toks[8]});
    const auto odds_ratio = std::stod(std::string{toks[9]});
    const auto omega = std::stod(std::string{toks[10]});

    constexpr auto max_coord = std::numeric_limits<std::uint32_t>::max();
    if (start1 > max_coord) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "start1 cannot be represented using 32 bit integers: value is too large: start1={}"),
          start1));
    }
    if (end1 > max_coord) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "end1 cannot be represented using 32 bit integers: value is too large: end1={}"),
          end1));
    }

    if (start2 > max_coord) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "start2 cannot be represented using 32 bit integers: value is too large: start1={}"),
          start2));
    }
    if (end2 > max_coord) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "end2 cannot be represented using 32 bit integers: value is too large: end1={}"),
          end2));
    }

    return {std::move(chrom1_ptr),
            static_cast<std::uint32_t>(start1),
            static_cast<std::uint32_t>(end1),
            std::move(chrom2_ptr),
            static_cast<std::uint32_t>(start2),
            static_cast<std::uint32_t>(end2),
            pvalue,
            pvalue,
            obs,
            exp,
            std::log2(odds_ratio) - std::log2(omega),
            odds_ratio,
            omega};
  }
};

int run_nchg_filter(const FilterConfig& c) {
  std::vector<Stats> records{};
  std::ifstream fs{};
  std::size_t i = 1;
  try {
    fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
    fs.open(c.path);

    ChromosomeSet chroms{};

    std::string buffer{};
    for (; std::getline(fs, buffer); ++i) {
      if (buffer.find('#') == 0) {
        continue;
      }
      records.emplace_back(Stats::parse(chroms, buffer));
    }
  } catch (const std::exception& e) {
    if (!fs.eof()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("an error occurred while parsing line {} from {}: {}"), i, c.path, e.what()));
    }
  }

  if (c.write_header) {
    print_header();
  }

  BH_FDR bh(std::move(records));
  for (const auto& s : bh.correct([](Stats& s) -> double& { return s.pvalue_corrected; })) {
    if (!c.drop_non_significant ||
        (s.pvalue_corrected <= c.fdr && std::abs(s.log_ratio) >= c.log_ratio)) {
      fmt::print(FMT_COMPILE("{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\t"
                             "{}\n"),
                 *s.chrom1, s.start1, s.end1, *s.chrom2, s.start2, s.end2, s.pvalue,
                 s.pvalue_corrected, s.observed_count, s.expected_count, s.log_ratio, s.odds_ratio,
                 s.omega);
    }
  }
  return 0;
}

}  // namespace nchg
