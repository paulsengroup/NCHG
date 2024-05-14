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

#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/io/concurrency.h>
#include <arrow/io/file.h>
#include <arrow/record_batch.h>
#include <arrow/util/thread_pool.h>
#include <blockingconcurrentqueue.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_writer.h>
#include <parquet/stream_reader.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <hictk/file.hpp>
#include <hictk/fmt/pixel.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <variant>

#include "nchg/common.hpp"
#include "nchg/config.hpp"
#include "nchg/io.hpp"
#include "nchg/k_merger.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {

[[nodiscard]] static std::pair<std::vector<ParquetStatsFile<NCHGResult>::iterator>,
                               std::vector<ParquetStatsFile<NCHGResult>::iterator>>
init_file_iterators(const std::filesystem::path &prefix, const hictk::Reference &chroms) {
  std::vector<ParquetStatsFile<NCHGResult>::iterator> heads{};
  std::vector<ParquetStatsFile<NCHGResult>::iterator> tails{};

  SPDLOG_INFO(FMT_STRING("enumerating chrom-chrom tables under prefix {}..."), prefix);
  for (std::uint32_t chrom1_id = 0; chrom1_id < chroms.size(); ++chrom1_id) {
    const auto &chrom1 = chroms.at(chrom1_id);
    if (chrom1.is_all()) {
      break;
    }
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chroms.size(); ++chrom2_id) {
      const auto &chrom2 = chroms.at(chrom2_id);
      const auto path = fmt::format(FMT_STRING("{}.{}.{}.parquet"), prefix.string(), chrom1.name(),
                                    chrom2.name());
      if (std::filesystem::exists(path)) {
        ParquetStatsFile<NCHGResult> f(chroms, path);
        auto first = f.begin();
        auto last = f.end();
        if (first != last) {
          heads.emplace_back(std::move(first));
          tails.emplace_back(std::move(last));
        }
      }
    }
  }

  SPDLOG_INFO(FMT_STRING("enumerated {} non-empty tables"), heads.size());

  return std::make_pair(std::move(heads), std::move(tails));
}

int run_nchg_merge(const MergeConfig &c) {
  const auto t0 = std::chrono::system_clock::now();

  const auto path_to_chrom_sizes =
      fmt::format(FMT_STRING("{}.chrom.sizes"), c.input_prefix.string());

  SPDLOG_INFO(FMT_STRING("reading chromosomes from \"{}\"..."), path_to_chrom_sizes);
  const auto chroms = hictk::Reference::from_chrom_sizes(path_to_chrom_sizes);
  SPDLOG_INFO(FMT_STRING("read {} chromosomes!"), chroms.size());

  moodycamel::BlockingConcurrentQueue<NCHGResult> queue(64 * 1024);

  BS::thread_pool tpool(static_cast<BS::concurrency_t>(std::min(2UL, c.threads)));

  std::atomic<bool> early_return{false};

  auto producer = tpool.submit_task([&]() {
    try {
      auto [heads, tails] = init_file_iterators(c.input_prefix, chroms);
      const KMerger merger(heads, tails);

      for (const auto &s : merger) {
        while (!queue.try_enqueue(s)) {
          if (early_return) {
            return;
          }
        }
      }

      NCHGResult s{};
      s.pval = -1;
      queue.enqueue(s);  // EOQ
    } catch (const std::exception &e) {
      early_return = true;
      throw std::runtime_error(
          fmt::format(FMT_STRING("an exception occurred in the producer thread: {}"), e.what()));
    } catch (...) {
      early_return = true;
      throw;
    }
  });

  auto consumer = tpool.submit_task([&]() {
    try {
      auto writer = init_parquet_file_writer(c.output_path, c.force, c.compression_method,
                                             c.compression_lvl, c.threads - 2);

      std::size_t records_processed = 0;
      NCHGResult buffer{};

      const std::size_t batch_size = 1'000'000;
      RecordBatchBuilder builder{};

      auto t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; true; ++i) {
        while (!queue.wait_dequeue_timed(buffer, std::chrono::milliseconds(10))) {
          if (early_return) {
            return records_processed;
          }
        }
        if (buffer.pval == -1) {  // EOQ
          break;
        }
        if (builder.size() == batch_size) {
          builder.write(*writer);
        }
        builder.append(buffer);
        ++records_processed;

        if (i == 10'000'000) {
          const auto t1 = std::chrono::steady_clock::now();
          const auto delta =
              static_cast<double>(
                  std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
              1000.0;

          SPDLOG_INFO(FMT_STRING("merging {:.0f} records/s..."), static_cast<double>(i) / delta);
          t0 = t1;
          i = 0;
        }
      }

      if (builder.size() != 0) {
        builder.write(*writer);
      }

      return records_processed;
    } catch (const std::exception &e) {
      early_return = true;
      throw std::runtime_error(
          fmt::format(FMT_STRING("an exception occurred in the consumer thread: {}"), e.what()));
    } catch (...) {
      early_return = true;
      throw;
    }
  });

  producer.get();
  const auto interactions_processed = consumer.get();

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  SPDLOG_INFO(FMT_STRING("Processed {} records in {}s!"), interactions_processed,
              static_cast<double>(delta) / 1000.0);

  return 0;
}

}  // namespace nchg
