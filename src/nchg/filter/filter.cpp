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
#include <nchg/nchg.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "nchg/fdr.hpp"
#include "nchg/io.hpp"
#include "nchg/tools.hpp"

namespace nchg {

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

struct CorrectedPvalue {
  std::size_t id{};

  double pvalue{};
  double pvalue_corrected{};
};

using ChromPair = std::pair<std::string, std::string>;

[[nodiscard]] static auto read_records(const std::filesystem::path& path) {
  SPDLOG_INFO(FMT_STRING("parsing records from {}"), path);
  phmap::btree_map<ChromPair, std::vector<CorrectedPvalue>> records{};

  ParquetStatsFile<NCHGResult> f(nullptr, path);

  std::size_t i = 0;
  for (const auto& record : f) {
    CorrectedPvalue pv{
        i++,
        record.pval,
        1.0,
    };

    auto [it, inserted] =
        records.try_emplace(std::make_pair(std::string{record.pixel.coords.bin1.chrom().name()},
                                           std::string{record.pixel.coords.bin2.chrom().name()}),
                            std::vector<CorrectedPvalue>{pv});

    if (!inserted) {
      it->second.emplace_back(pv);
    }
  }

  return records;
}

[[nodiscard]] static std::vector<CorrectedPvalue> correct_pvalues_chrom_chrom(
    const phmap::btree_map<ChromPair, std::vector<CorrectedPvalue>>& records) {
  SPDLOG_INFO(
      FMT_STRING("proceeding to correct pvalues for individual chromosme pairs separately..."));
  std::vector<CorrectedPvalue> corrected_records{};
  BH_FDR<CorrectedPvalue> bh{};
  for (const auto& [cp, values] : records) {
    SPDLOG_INFO(FMT_STRING("processing {}:{} values..."), cp.first, cp.second);
    bh.clear();
    bh.add_records(values.begin(), values.end());
    auto chrom_chrom_corrected_records =
        bh.correct([](CorrectedPvalue& s) -> double& { return s.pvalue_corrected; });
    corrected_records.insert(corrected_records.end(),
                             std::make_move_iterator(chrom_chrom_corrected_records.begin()),
                             std::make_move_iterator(chrom_chrom_corrected_records.end()));
  }
  return corrected_records;
}

[[nodiscard]] static std::vector<CorrectedPvalue> correct_pvalues_cis_trans(
    const phmap::btree_map<ChromPair, std::vector<CorrectedPvalue>>& records) {
  SPDLOG_INFO(FMT_STRING("proceeding to correct pvalues for cis and trans matrices separately..."));
  BH_FDR<CorrectedPvalue> bh_cis{};
  BH_FDR<CorrectedPvalue> bh_trans{};
  for (const auto& [chroms, values] : records) {
    if (chroms.first == chroms.second) {
      bh_cis.add_records(values.begin(), values.end());
    } else {
      bh_trans.add_records(values.begin(), values.end());
    }
  }

  auto corrected_records =
      bh_cis.correct([](CorrectedPvalue& s) -> double& { return s.pvalue_corrected; });
  auto corrected_records_trans =
      bh_trans.correct([](CorrectedPvalue& s) -> double& { return s.pvalue_corrected; });
  corrected_records.insert(corrected_records.end(),
                           std::make_move_iterator(corrected_records_trans.begin()),
                           std::make_move_iterator(corrected_records_trans.end()));
  return corrected_records;
}

[[nodiscard]] static std::vector<CorrectedPvalue> correct_pvalues(
    const phmap::btree_map<ChromPair, std::vector<CorrectedPvalue>>& records,
    const FilterConfig& c) {
  if (c.correct_chrom_chrom_separately) {
    return correct_pvalues_chrom_chrom(records);
  }

  if (c.correct_cis_trans_separately) {
    return correct_pvalues_cis_trans(records);
  }

  SPDLOG_INFO(FMT_STRING("proceeding to correct all pvalues in one go"));
  BH_FDR<CorrectedPvalue> bh{};
  for (const auto& [_, values] : records) {
    bh.add_records(values.begin(), values.end());
  }

  return bh.correct([](CorrectedPvalue& s) -> double& { return s.pvalue_corrected; });
}

struct NCHGFilterResult {
  hictk::Pixel<std::uint32_t> pixel{};
  double expected{};
  double pval{};
  double pval_corrected{};
  double log_ratio{};
  double odds_ratio{};
  double omega{};
};

int run_nchg_filter(const FilterConfig& c) {
  auto corrected_pvalues = correct_pvalues(read_records(c.input_path), c);

  std::sort(corrected_pvalues.begin(), corrected_pvalues.end(),
            [&](const auto& r1, const auto& r2) { return r1.id < r2.id; });

  auto writer = init_parquet_file_writer(c.output_path, c.force, c.compression_method,
                                         c.compression_lvl, c.threads - 1);

  const std::size_t batch_size = 1'000'000;
  RecordBatchBuilder builder{};

  std::size_t i = 0;
  for (const auto& r : ParquetStatsFile<NCHGResult>(nullptr, c.input_path)) {
    const auto pvalue_corrected = corrected_pvalues[i++].pvalue_corrected;
    const auto log_ratio = std::log2(r.odds_ratio) - std::log2(r.omega);

    if (!c.drop_non_significant || (pvalue_corrected <= c.fdr && log_ratio >= c.log_ratio)) {
      NCHGFilterResult res{r.pixel,   r.expected,   r.pval, pvalue_corrected,
                           log_ratio, r.odds_ratio, r.omega};
      if (builder.size() == batch_size) {
        builder.write(*writer);
      }
      builder.append(res);
    }

    if (builder.size() != 0) {
      builder.write(*writer);
    }
  }
  return 0;
}

}  // namespace nchg
