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
#include <readerwriterqueue/readerwriterqueue.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
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

using ChromPair = std::pair<std::string, std::string>;

[[nodiscard]] static auto read_records(const std::filesystem::path& path) {
  SPDLOG_INFO(FMT_STRING("reading records from {}"), path);
  phmap::btree_map<ChromPair, std::vector<double>> records{};

  ParquetStatsFile<NCHGResult> f(nullptr, path);

  for (const auto& record : f) {
    auto [it, inserted] =
        records.try_emplace(std::make_pair(std::string{record.pixel.coords.bin1.chrom().name()},
                                           std::string{record.pixel.coords.bin2.chrom().name()}),
                            std::vector<double>{record.pval});

    if (!inserted) {
      it->second.emplace_back(record.pval);
    }
  }

  return records;
}

[[nodiscard]] static std::vector<double> correct_pvalues_chrom_chrom(
    const phmap::btree_map<ChromPair, std::vector<double>>& records) {
  SPDLOG_INFO(
      FMT_STRING("proceeding to correct pvalues for individual chromosme pairs separately..."));
  std::vector<double> corrected_records{};
  BH_FDR<double> bh{};
  for (const auto& [cp, values] : records) {
    SPDLOG_INFO(FMT_STRING("processing {}:{} values..."), cp.first, cp.second);
    bh.clear();
    bh.add_records(values.begin(), values.end());
    auto chrom_chrom_corrected_records = bh.correct();
    corrected_records.insert(corrected_records.end(),
                             std::make_move_iterator(chrom_chrom_corrected_records.begin()),
                             std::make_move_iterator(chrom_chrom_corrected_records.end()));
  }
  return corrected_records;
}

[[nodiscard]] static std::vector<double> correct_pvalues_cis_trans(
    const phmap::btree_map<ChromPair, std::vector<double>>& records) {
  SPDLOG_INFO(FMT_STRING("proceeding to correct pvalues for cis and trans matrices separately..."));
  BH_FDR<double> bh_cis{};
  BH_FDR<double> bh_trans{};
  for (const auto& [chroms, values] : records) {
    if (chroms.first == chroms.second) {
      bh_cis.add_records(values.begin(), values.end());
    } else {
      bh_trans.add_records(values.begin(), values.end());
    }
  }

  auto corrected_records = bh_cis.correct();
  auto corrected_records_trans = bh_trans.correct();
  corrected_records.insert(corrected_records.end(),
                           std::make_move_iterator(corrected_records_trans.begin()),
                           std::make_move_iterator(corrected_records_trans.end()));
  return corrected_records;
}

[[nodiscard]] static std::vector<double> correct_pvalues(
    const phmap::btree_map<ChromPair, std::vector<double>>& records, const FilterConfig& c) {
  if (c.correct_chrom_chrom_separately) {
    return correct_pvalues_chrom_chrom(records);
  }

  if (c.correct_cis_trans_separately) {
    return correct_pvalues_cis_trans(records);
  }

  SPDLOG_INFO(FMT_STRING("proceeding to correct all pvalues in one go"));
  BH_FDR<double> bh{};
  for (const auto& [_, values] : records) {
    bh.add_records(values.begin(), values.end());
  }

  return bh.correct();
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

static_assert(has_pval_corrected<NCHGFilterResult>::value);
static_assert(has_log_ratio<NCHGFilterResult>::value);

int run_nchg_filter(const FilterConfig& c) {
  auto corrected_pvalues = correct_pvalues(read_records(c.input_path), c);

  BS::thread_pool tpool(2);

  moodycamel::BlockingReaderWriterQueue<NCHGFilterResult> queue{64 * 1024};
  std::atomic<bool> early_return{};

  auto producer = tpool.submit_task([&]() {
    std::size_t i = 0;
    try {
      for (const auto& r : ParquetStatsFile<NCHGResult>(nullptr, c.input_path)) {
        if (early_return) {
          return;
        }

        const auto pvalue_corrected = corrected_pvalues[i++];

        if (!c.drop_non_significant || (pvalue_corrected <= c.fdr && r.log_ratio >= c.log_ratio)) {
          const NCHGFilterResult res{r.pixel,     r.expected,   r.pval, pvalue_corrected,
                                     r.log_ratio, r.odds_ratio, r.omega};
          while (!queue.try_enqueue(res)) {
            if (early_return) {
              return;
            }
          }
        }
      }

      NCHGFilterResult eoq{};
      eoq.pval = -1;
      queue.enqueue(eoq);

    } catch (const std::exception& e) {
      early_return = true;
      throw std::runtime_error(
          fmt::format(FMT_STRING("an exception occurred in producer thread: {}"), e.what()));
    } catch (...) {
      early_return = true;
      throw;
    }
  });

  auto consumer = tpool.submit_task([&]() {
    SPDLOG_INFO(FMT_STRING("writing records to output file {}"), c.output_path);
    auto writer = init_parquet_file_writer<NCHGFilterResult>(
        c.output_path, c.force, c.compression_method, c.compression_lvl, c.threads - 2);

    const std::size_t batch_size = 1'000'000;
    RecordBatchBuilder builder{};

    NCHGFilterResult res{};

    while (!early_return) {
      while (!queue.wait_dequeue_timed(res, std::chrono::milliseconds(20))) {
        if (early_return) {
          return;
        }
      }

      if (res.pval == -1) {
        // EOQ
        if (builder.size() != 0) {
          builder.write(*writer);
        }
        return;
      }

      if (builder.size() == batch_size) {
        builder.write(*writer);
        builder.reset();
      }
      builder.append(res);
    }
  });

  producer.wait();
  consumer.wait();

  return 0;
}

}  // namespace nchg
