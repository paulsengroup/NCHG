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

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <hictk/chromosome.hpp>
#include <hictk/numeric_utils.hpp>
#include <iostream>
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

template <typename N>
[[nodiscard]] static N parse_numeric(std::string_view tok) {
  return hictk::internal::parse_numeric_or_throw<N>(tok);
}

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

  [[nodiscard]] static Stats parse(std::vector<std::string_view>& toks, ChromosomeSet& chromosomes,
                                   std::string_view s, char sep = '\t') {
    if (s.empty()) {
      throw std::runtime_error("found an empty line");
    }

    if (s.back() == '\r') {
      s = s.substr(0, s.size() - 1);
    }

    toks.clear();
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
    const auto start1 = parse_numeric<std::uint64_t>(toks[1]);
    const auto end1 = parse_numeric<std::uint64_t>(toks[2]);

    if (chromosomes.find(toks[3]) == chromosomes.end()) {
      chromosomes.emplace(std::make_shared<std::string>(std::string{toks[3]}));
    }
    auto chrom2_ptr = *chromosomes.find(toks[3]);
    const auto start2 = parse_numeric<std::uint64_t>(toks[4]);
    const auto end2 = parse_numeric<std::uint64_t>(toks[5]);

    const auto pvalue = parse_numeric<double>(toks[6]);
    const auto obs = parse_numeric<std::uint64_t>(toks[7]);
    const auto exp = parse_numeric<double>(toks[8]);
    const auto odds_ratio = parse_numeric<double>(toks[9]);
    const auto omega = parse_numeric<double>(toks[10]);

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

    const auto log_ratio = std::log2(odds_ratio) - std::log2(omega);

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
            std::isfinite(log_ratio) ? log_ratio : 0.0,
            odds_ratio,
            omega};
  }

  [[nodiscard]] bool operator<(const Stats& other) const noexcept {
    if (chrom1 != other.chrom1) {
      return *chrom1 < *other.chrom1;
    }
    if (start1 != other.start1) {
      return start1 < other.start1;
    }
    if (chrom2 != other.chrom2) {
      return *chrom2 < *other.chrom2;
    }
    return start2 < other.start2;
  }
};

using ChromPair = std::pair<std::string, std::string>;

[[nodiscard]] std::istream* open_file(const std::filesystem::path& path, std::ifstream& ifs) {
  if (path.empty() || path == "-") {
    return &std::cin;
  }

  ifs.exceptions(ifs.exceptions() | std::ios::failbit | std::ios::badbit);
  ifs.open(path);

  return &ifs;
}

[[nodiscard]] static phmap::btree_map<ChromPair, std::vector<Stats>> parse_records(
    const std::filesystem::path& path) {
  SPDLOG_INFO(FMT_STRING("parsing records from {}"), path.empty() || path == "-" ? "stdin" : path);
  phmap::btree_map<ChromPair, std::vector<Stats>> records{};

  std::ifstream ifs{};
  std::istream* is = &std::cin;
  std::size_t i = 1;
  try {
    ChromosomeSet chroms{};

    if (!path.empty() && path != "-") {
      is = open_file(path, ifs);
    }

    std::vector<std::string_view> tok_buffer{};
    std::string buffer{};
    for (; std::getline(*is, buffer); ++i) {
      if (i == 1 && buffer.find("chrom1") == 0) {
        continue;
      }
      auto record = Stats::parse(tok_buffer, chroms, buffer);
      auto [it, inserted] = records.try_emplace(std::make_pair(*record.chrom1, *record.chrom2),
                                                std::vector<Stats>{record});
      if (!inserted) {
        it->second.emplace_back(std::move(record));
      }
    }
  } catch (const std::exception& e) {
    if (!is->eof()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("an error occurred while parsing line {} from {}: {}"), i,
                      is == &std::cin ? "stdin" : path, e.what()));
    }
  }

  return records;
}

[[nodiscard]] static std::vector<Stats> correct_pvalues_chrom_chrom(
    const phmap::btree_map<ChromPair, std::vector<Stats>>& records) {
  SPDLOG_INFO(
      FMT_STRING("proceeding to correct pvalues for individual chromosme pairs separately..."));
  std::vector<Stats> corrected_records{};
  BH_FDR<Stats> bh{};
  for (const auto& [cp, values] : records) {
    SPDLOG_INFO(FMT_STRING("processing {}:{} values..."), cp.first, cp.second);
    bh.clear();
    bh.add_records(values.begin(), values.end());
    auto chrom_chrom_corrected_records =
        bh.correct([](Stats& s) -> double& { return s.pvalue_corrected; });
    corrected_records.insert(corrected_records.end(),
                             std::make_move_iterator(chrom_chrom_corrected_records.begin()),
                             std::make_move_iterator(chrom_chrom_corrected_records.end()));
  }
  return corrected_records;
}

[[nodiscard]] static std::vector<Stats> correct_pvalues_cis_trans(
    const phmap::btree_map<ChromPair, std::vector<Stats>>& records) {
  SPDLOG_INFO(FMT_STRING("proceeding to correct pvalues for cis and trans matrices separately..."));
  BH_FDR<Stats> bh_cis{};
  BH_FDR<Stats> bh_trans{};
  for (const auto& [chroms, values] : records) {
    if (chroms.first == chroms.second) {
      bh_cis.add_records(values.begin(), values.end());
    } else {
      bh_trans.add_records(values.begin(), values.end());
    }
  }

  auto corrected_records = bh_cis.correct([](Stats& s) -> double& { return s.pvalue_corrected; });
  auto corrected_records_trans =
      bh_trans.correct([](Stats& s) -> double& { return s.pvalue_corrected; });
  corrected_records.insert(corrected_records.end(),
                           std::make_move_iterator(corrected_records_trans.begin()),
                           std::make_move_iterator(corrected_records_trans.end()));
  return corrected_records;
}

[[nodiscard]] static std::vector<Stats> correct_pvalues(
    const phmap::btree_map<ChromPair, std::vector<Stats>>& records, const FilterConfig& c) {
  if (c.correct_chrom_chrom_separately) {
    return correct_pvalues_chrom_chrom(records);
  }

  if (c.correct_cis_trans_separately) {
    return correct_pvalues_cis_trans(records);
  }

  SPDLOG_INFO(FMT_STRING("proceeding to correct all pvalues in one go"));
  BH_FDR<Stats> bh{};
  for (const auto& [_, values] : records) {
    bh.add_records(values.begin(), values.end());
  }

  return bh.correct([](Stats& s) -> double& { return s.pvalue_corrected; });
}

int run_nchg_filter(const FilterConfig& c) {
  auto records = correct_pvalues(parse_records(c.path), c);

  if (c.sorted) {
    SPDLOG_INFO("sorting records...");
    std::sort(records.begin(), records.end());
  }

  if (c.write_header) {
    SPDLOG_INFO("printing the file header...");
    print_header();
  }

  for (const auto& s : records) {
    if (!c.drop_non_significant || (s.pvalue_corrected <= c.fdr && s.log_ratio >= c.log_ratio)) {
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
