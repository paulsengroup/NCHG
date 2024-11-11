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

// clang-format off
#include <queue>


#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <hictk/reference.hpp>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
#include <moodycamel/blockingconcurrentqueue.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <future>
#include <stdexcept>
#include <utility>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/k_merger.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/io.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

[[nodiscard]] static std::pair<std::vector<ParquetStatsFile::iterator<NCHGResult>>,
                               std::vector<ParquetStatsFile::iterator<NCHGResult>>>
init_file_iterators(const std::filesystem::path &prefix, const hictk::Reference &chroms) {
  std::vector<ParquetStatsFile::iterator<NCHGResult>> heads{};
  std::vector<ParquetStatsFile::iterator<NCHGResult>> tails{};

  SPDLOG_INFO("enumerating chrom-chrom tables under prefix {}...", prefix);
  for (std::uint32_t chrom1_id = 0; chrom1_id < chroms.size(); ++chrom1_id) {
    const auto &chrom1 = chroms.at(chrom1_id);
    if (chrom1.is_all()) [[unlikely]] {
      break;
    }
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chroms.size(); ++chrom2_id) {
      const auto &chrom2 = chroms.at(chrom2_id);
      const auto path =
          fmt::format("{}.{}.{}.parquet", prefix.string(), chrom1.name(), chrom2.name());
      if (std::filesystem::exists(path)) {
        ParquetStatsFile f(path, ParquetStatsFile::RecordType::NCHGCompute);
        auto first = f.begin<NCHGResult>();
        auto last = f.end<NCHGResult>();
        if (first != last) [[likely]] {
          heads.emplace_back(std::move(first));
          tails.emplace_back(std::move(last));
        }
      }
    }
  }

  if (heads.empty()) {
    throw std::runtime_error(fmt::format("unable to find any table under prefix {}", prefix));
  }

  SPDLOG_INFO("enumerated {} non-empty tables", heads.size());

  return {heads, tails};
}

using RecordQueue = moodycamel::BlockingConcurrentQueue<NCHGResult>;

[[nodiscard]] static std::size_t producer_fx(const hictk::Reference &chromosomes,
                                             const std::filesystem::path &input_prefix,
                                             RecordQueue &queue, std::atomic<bool> &early_return) {
  try {
    std::size_t records_enqueued{};
    auto [heads, tails] = init_file_iterators(input_prefix, chromosomes);
    const KMerger merger(heads, tails);

    for (const auto &s : merger) {
      while (!queue.try_enqueue(s)) [[unlikely]] {
        if (early_return) [[unlikely]] {
          return records_enqueued;
        }
      }
      ++records_enqueued;
    }

    NCHGResult s{};
    s.pval = -1;
    queue.enqueue(s);  // EOQ

    return records_enqueued;
  } catch (const std::exception &e) {
    early_return = true;
    throw std::runtime_error(fmt::format("an exception occurred in producer thread: {}", e.what()));
  } catch (...) {
    SPDLOG_ERROR("an unknown exception occurred in producer thread");
    early_return = true;
    throw;
  }
}

[[nodiscard]] static std::size_t consumer_fx(const MergeConfig &c,
                                             const hictk::Reference &chromosomes,
                                             RecordQueue &queue, std::atomic<bool> &early_return) {
  try {
    auto writer = init_parquet_file_writer<NCHGResult>(chromosomes, c.output_path, c.force,
                                                       c.compression_method, c.compression_lvl,
                                                       c.threads - 2);

    std::size_t records_dequeued = 0;
    NCHGResult buffer{};

    const std::size_t batch_size = 1'000'000;
    RecordBatchBuilder builder(chromosomes);

    auto t1 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; true; ++i) {
      while (!queue.wait_dequeue_timed(buffer, std::chrono::milliseconds(10))) [[unlikely]] {
        if (early_return) [[unlikely]] {
          return records_dequeued;
        }
      }
      if (buffer.pval == -1) [[unlikely]] {  // EOQ
        break;
      }

      if (builder.size() == batch_size) [[unlikely]] {
        builder.write(*writer);
      }
      builder.append(buffer);
      ++records_dequeued;

      if (i == 10'000'000) [[unlikely]] {
        const auto t2 = std::chrono::steady_clock::now();
        const auto delta =
            static_cast<double>(
                std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()) /
            1000.0;

        SPDLOG_INFO("merging {:.0f} records/s...", static_cast<double>(i) / delta);
        t1 = t2;
        i = 0;
      }
    }

    if (builder.size() != 0) {
      builder.write(*writer);
    }

    return records_dequeued;
  } catch (const std::exception &e) {
    early_return = true;
    throw std::runtime_error(fmt::format("an exception occurred in consumer thread: {}", e.what()));
  } catch (...) {
    SPDLOG_ERROR("an unknown exception occurred in consumer thread");
    early_return = true;
    throw;
  }
}

int run_command(const MergeConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();

  const auto path_to_chrom_sizes = fmt::format("{}.chrom.sizes", c.input_prefix.string());

  SPDLOG_INFO("reading chromosomes from \"{}\"...", path_to_chrom_sizes);
  const auto chroms = hictk::Reference::from_chrom_sizes(path_to_chrom_sizes);
  SPDLOG_INFO("read {} chromosomes!", chroms.size());

  moodycamel::BlockingConcurrentQueue<NCHGResult> queue(64 * 1024);
  std::atomic<bool> early_return{false};

  auto producer = std::async(std::launch::deferred, [&] {
    SPDLOG_DEBUG("spawning producer thread...");
    return producer_fx(chroms, c.input_prefix, queue, early_return);
  });

  auto consumer = std::async(std::launch::async, [&] {
    SPDLOG_DEBUG("spawning consumer thread...");
    return consumer_fx(c, chroms, queue, early_return);
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
