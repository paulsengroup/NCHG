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
#include <parallel_hashmap/btree.h>
#include <readerwriterqueue/readerwriterqueue.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <future>
#include <hictk/chromosome.hpp>
#include <hictk/numeric_utils.hpp>
#include <memory>
#include <nchg/nchg.hpp>
#include <numeric>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "nchg/fdr.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/io.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

struct SharedPtrStringCmp {
  using is_transparent = void;
  constexpr bool operator()(const std::shared_ptr<std::string>& a,
                            const std::shared_ptr<std::string>& b) const noexcept {
    assert(!!a);
    assert(!!b);
    return *a < *b;
  }

  constexpr bool operator()(const std::string& a,
                            const std::shared_ptr<std::string>& b) const noexcept {
    assert(!!b);
    return a < *b;
  }

  constexpr bool operator()(const std::shared_ptr<std::string>& a,
                            const std::string& b) const noexcept {
    assert(!!a);
    return *a < b;
  }

  constexpr bool operator()(std::string_view a,
                            const std::shared_ptr<std::string>& b) const noexcept {
    assert(!!b);
    return a < *b;
  }

  constexpr bool operator()(const std::shared_ptr<std::string>& a,
                            std::string_view b) const noexcept {
    assert(!!a);
    return *a < b;
  }
};

struct StringPairCmp {
  using is_transparent = void;
  constexpr bool operator()(const std::pair<std::string, std::string>& a,
                            const std::pair<std::string, std::string>& b) const noexcept {
    return a < b;
  }

  constexpr bool operator()(const std::pair<std::string_view, std::string_view>& a,
                            const std::pair<std::string_view, std::string_view>& b) const noexcept {
    return a < b;
  }

  constexpr bool operator()(const std::pair<std::string, std::string>& a,
                            const std::pair<std::string_view, std::string_view>& b) const noexcept {
    if (a.first == b.first) {
      return a.second < b.second;
    }
    return a.first < b.first;
  }

  constexpr bool operator()(const std::pair<std::string_view, std::string_view>& a,
                            const std::pair<std::string, std::string>& b) const noexcept {
    if (a.first == b.first) {
      return a.second < b.second;
    }
    return a.first < b.first;
  }
};

template <typename N>
[[nodiscard]] static N parse_numeric(std::string_view tok) {
  return hictk::internal::parse_numeric_or_throw<N>(tok);
}

using ChromPair = std::pair<std::string, std::string>;

struct PValue {
  std::size_t i{};
  double pvalue{};
};

using ChromChromPvalueMap = phmap::btree_map<ChromPair, std::vector<PValue>, StringPairCmp>;
using FlatPvalueMap = phmap::flat_hash_map<std::size_t, double>;

[[nodiscard]] static ChromChromPvalueMap read_records(const std::filesystem::path& path) {
  SPDLOG_INFO("reading records from {}", path);
  ChromChromPvalueMap records{};

  std::size_t i = 0;
  ParquetStatsFile f{path, ParquetStatsFile::RecordType::NCHGCompute};
  std::for_each(f.begin<NCHGResult>(), f.end<NCHGResult>(), [&](const NCHGResult& record) {
    const auto& chrom1 = record.pixel.coords.bin1.chrom();
    const auto& chrom2 = record.pixel.coords.bin2.chrom();

    auto it = records.find(std::make_pair(chrom1.name(), chrom2.name()));
    if (it == records.end()) [[unlikely]] {
      it = records
               .emplace(std::make_pair(std::string{chrom1.name()}, std::string{chrom2.name()}),
                        std::vector<PValue>{})
               .first;
    }

#if defined(__apple_build_version__) && __apple_build_version__ < 16000000
    it->second.emplace_back(PValue{i, record.pval});
#else
    it->second.emplace_back(i, record.pval);
#endif
    ++i;
  });

  return records;
}

[[nodiscard]] static FlatPvalueMap alloc_pvalue_hashmap(const ChromChromPvalueMap& records) {
  auto sizes =
      records | std::views::values | std::views::transform([](const auto& v) { return v.size(); });
  const auto num_records = std::accumulate(sizes.begin(), sizes.end(), 0uz);

  return FlatPvalueMap{num_records};
}

[[nodiscard]] static FlatPvalueMap correct_pvalues_chrom_chrom(const ChromChromPvalueMap& records) {
  SPDLOG_INFO("proceeding to correct pvalues for individual chromosome pairs separately...");

  const auto t0 = std::chrono::steady_clock::now();

  auto corrected_records = alloc_pvalue_hashmap(records);

  BH_FDR<PValue> bh;
  for (const auto& [cp, values] : records) {
    SPDLOG_INFO("processing {}:{} values...", cp.first, cp.second);
    const auto t1 = std::chrono::steady_clock::now();
    bh.clear();
    bh.add_records(values);

    for (const auto& record : bh.correct([](auto& r) -> double& { return r.pvalue; })) {
      corrected_records.emplace(record.i, record.pvalue);
    }
    const auto t2 = std::chrono::steady_clock::now();
    SPDLOG_INFO("correcting pvalues for {}:{} took {}", cp.first, cp.second,
                format_duration(t2 - t1));
  }

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("corrected {} pvalues in {}", corrected_records.size(), format_duration(t1 - t0));

  return corrected_records;
}

[[nodiscard]] static FlatPvalueMap correct_pvalues_cis_trans(const ChromChromPvalueMap& records) {
  SPDLOG_INFO("proceeding to correct pvalues for cis and trans matrices separately...");

  const auto t0 = std::chrono::steady_clock::now();

  auto corrected_records = alloc_pvalue_hashmap(records);

  BH_FDR<PValue> bh_cis{};
  BH_FDR<PValue> bh_trans{};
  for (const auto& [chroms, values] : records) {
    if (chroms.first == chroms.second) {
      bh_cis.add_records(values);
    } else {
      bh_trans.add_records(values);
    }
  }

  for (const auto& record : bh_cis.correct([](auto& r) -> double& { return r.pvalue; })) {
    corrected_records.emplace(record.i, record.pvalue);
  }

  for (const auto& record : bh_trans.correct([](auto& r) -> double& { return r.pvalue; })) {
    corrected_records.emplace(record.i, record.pvalue);
  }

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("corrected {} pvalues in {}", corrected_records.size(), format_duration(t1 - t0));

  return corrected_records;
}

[[nodiscard]] static FlatPvalueMap correct_pvalues(const ChromChromPvalueMap& records,
                                                   const FilterConfig& c) {
  if (c.correct_chrom_chrom_separately) {
    return correct_pvalues_chrom_chrom(records);
  }

  if (c.correct_cis_trans_separately) {
    return correct_pvalues_cis_trans(records);
  }

  SPDLOG_INFO("proceeding to correct all pvalues in one go");
  const auto t0 = std::chrono::steady_clock::now();

  BH_FDR<PValue> bh{};
  for (const auto& values : records | std::views::values) {
    bh.add_records(values);
  }

  auto corrected_records = alloc_pvalue_hashmap(records);

  for (const auto& record : bh.correct([](auto& r) -> double& { return r.pvalue; })) {
    corrected_records.emplace(std::make_pair(record.i, record.pvalue));
  }

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("corrected {} pvalues in {}", corrected_records.size(), format_duration(t1 - t0));
  return corrected_records;
}

struct NCHGFilterResult {
  hictk::Pixel<std::uint64_t> pixel{};
  double expected{};
  double pval{};
  double pval_corrected{};
  double log_ratio{};
  double odds_ratio{};
  double omega{};

  NCHGFilterResult() = default;
  NCHGFilterResult(const NCHGResult& r, double pval_corrected_)
      : pixel(r.pixel),
        expected(r.expected),
        pval(r.pval),
        pval_corrected(pval_corrected_),
        log_ratio(r.log_ratio),
        odds_ratio(r.odds_ratio),
        omega(r.omega) {}
};

static_assert(has_pval_corrected<NCHGFilterResult>::value);
static_assert(has_log_ratio<NCHGFilterResult>::value);

using RecordQueue = moodycamel::BlockingReaderWriterQueue<NCHGFilterResult>;

[[nodiscard]] static std::size_t producer_fx(const FilterConfig& c,
                                             const FlatPvalueMap& corrected_pvalues,
                                             RecordQueue& queue, std::atomic<bool>& early_return) {
  std::size_t records_enqueued{};
  try {
    ParquetStatsFile f{c.input_path, ParquetStatsFile::RecordType::NCHGCompute};
    auto first = f.begin<NCHGResult>();
    const auto last = f.end<NCHGResult>();
    for (; first != last; ++first) {
      if (early_return) [[unlikely]] {
        return records_enqueued;
      }

      const auto pvalue_corrected = corrected_pvalues.at(records_enqueued++);
      assert(first->pval <= pvalue_corrected);

      if (!c.drop_non_significant || (pvalue_corrected <= c.fdr && first->log_ratio >= c.log_ratio))
          [[likely]] {
        const NCHGFilterResult res{*first, pvalue_corrected};
        while (!queue.try_enqueue(res)) [[unlikely]] {
          if (early_return) [[unlikely]] {
            return records_enqueued;
          }
        }
      }
    }

    NCHGFilterResult eoq{};
    eoq.pval = -1;
    queue.enqueue(eoq);

    return records_enqueued;

  } catch (const std::exception& e) {
    early_return = true;
    throw std::runtime_error(fmt::format("an exception occurred in producer thread: {}", e.what()));
  } catch (...) {
    SPDLOG_ERROR("an unknown exception occurred in producer thread");
    early_return = true;
    throw;
  }
}

[[nodiscard]] static std::size_t consumer_fx(const FilterConfig& c, RecordQueue& queue,
                                             std::atomic<bool>& early_return) {
  SPDLOG_INFO("writing records to output file {}", c.output_path);
  std::size_t records_dequeued{};
  try {
    const auto chroms =
        *ParquetStatsFile(c.input_path, ParquetStatsFile::RecordType::NCHGCompute).chromosomes();

    const auto writer = init_parquet_file_writer<NCHGFilterResult>(
        chroms, c.output_path, c.force, c.compression_method, c.compression_lvl, c.threads - 2);

    constexpr std::size_t batch_size = 1'000'000;
    RecordBatchBuilder builder{chroms};

    NCHGFilterResult res{};

    while (!early_return) [[likely]] {
      while (!queue.wait_dequeue_timed(res, std::chrono::milliseconds(20))) [[unlikely]] {
        if (early_return) [[likely]] {
          return records_dequeued;
        }
      }

      if (res.pval == -1) [[unlikely]] {
        // EOQ
        if (builder.size() != 0) {
          builder.write(*writer);
        }
        return records_dequeued;
      }

      if (builder.size() == batch_size) [[unlikely]] {
        builder.write(*writer);
        builder.reset();
      }
      builder.append(res);
      ++records_dequeued;
    }

    return records_dequeued;
  } catch (const std::exception& e) {
    early_return = true;
    throw std::runtime_error(fmt::format("an exception occurred in consumer thread: {}", e.what()));
  } catch (...) {
    SPDLOG_ERROR("an unknown exception occurred in consumer thread");
    early_return = true;
    throw;
  }
}

int run_nchg_filter(const FilterConfig& c) {
  const auto t0 = std::chrono::steady_clock::now();
  const auto corrected_pvalues = correct_pvalues(read_records(c.input_path), c);

  RecordQueue queue{64uz * 1024uz};
  std::atomic<bool> early_return{};

  auto producer = std::async(std::launch::deferred, [&] {
    SPDLOG_DEBUG("spawning producer thread...");
    return producer_fx(c, corrected_pvalues, queue, early_return);
  });
  auto consumer = std::async(std::launch::async, [&] {
    SPDLOG_DEBUG("spawning consumer thread...");
    return consumer_fx(c, queue, early_return);
  });
  const auto records_enqueued = producer.get();
  const auto records_dequeued = consumer.get();

  if (records_enqueued != records_dequeued) {
    throw std::runtime_error(
        fmt::format("record queue is corrupted: not all records have been dequeued: "
                    "expected {} records, found {}",
                    records_enqueued, records_dequeued));
  }

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("processed {} records in {}!", records_dequeued, format_duration(t1 - t0));

  return 0;
}

}  // namespace nchg
