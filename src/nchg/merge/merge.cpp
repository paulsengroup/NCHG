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
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <future>
#include <hictk/reference.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "nchg/file_metadata.hpp"
#include "nchg/k_merger.hpp"
#include "nchg/nchg.hpp"
#include "nchg/parquet_stats_file_reader.hpp"
#include "nchg/parquet_stats_file_writer.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

using FileIteratorVect = std::vector<ParquetStatsFileReader::iterator<NCHGResult>>;
struct FileIteratorPairs {
  FileIteratorVect heads;
  FileIteratorVect tails;
};

[[nodiscard]] static std::string generate_chrom_sizes_name(
    const std::filesystem::path &input_prefix) {
  return fmt::format("{}.chrom.sizes", input_prefix.string());
}

[[nodiscard]] static std::string generate_report_name(const std::filesystem::path &input_prefix) {
  return fmt::format("{}.json", input_prefix);
}

[[nodiscard]] static FileIteratorPairs init_file_iterator_pairs(
    const NCHGResultMetadata &metadata, const std::filesystem::path &input_prefix) {
  try {
    assert(!metadata.records().empty());

    SPDLOG_INFO("enumerating tables based on file(s) listed in report file \"{}\"...",
                generate_report_name(input_prefix));

    const auto root_folder = input_prefix.has_parent_path() ? input_prefix.parent_path()
                                                            : std::filesystem::current_path();

    FileIteratorPairs result;
    auto &[heads, tails] = result;
    heads.reserve(metadata.records().size() - 1);
    tails.reserve(metadata.records().size() - 1);

    for (const auto &record : metadata.records()) {
      if (record.name.extension() != ".parquet") {
        continue;
      }

      try {
        ParquetStatsFileReader f(root_folder / record.name,
                                 ParquetStatsFileReader::RecordType::NCHGCompute);
        auto first = f.begin<NCHGResult>();
        auto last = f.end<NCHGResult>();
        if (first != last) [[likely]] {
          heads.emplace_back(std::move(first));
          tails.emplace_back(std::move(last));
        }
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format("failed to initialize iterators for file \"{}\": {}",
                                             root_folder / record.name, e.what()));
      }
    }

    return result;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to enumerate tables under prefix \"{}\": {}", input_prefix, e.what()));
  }
}

[[nodiscard]] static FileIteratorPairs init_file_iterator_pairs(
    const std::filesystem::path &input_prefix, const hictk::Reference &chroms) {
  SPDLOG_INFO("enumerating table based on chromosome(s) listed in \"{}\"...",
              generate_chrom_sizes_name(input_prefix));
  FileIteratorPairs result;
  auto &[heads, tails] = result;

  for (std::uint32_t chrom1_id = 0; chrom1_id < chroms.size(); ++chrom1_id) {
    const auto &chrom1 = chroms.at(chrom1_id);
    if (chrom1.is_all()) [[unlikely]] {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chroms.size(); ++chrom2_id) {
      const auto &chrom2 = chroms.at(chrom2_id);
      if (chrom2.is_all()) [[unlikely]] {
        continue;
      }

      const auto path =
          fmt::format("{}.{}.{}.parquet", input_prefix.string(), chrom1.name(), chrom2.name());
      if (std::filesystem::exists(path)) {
        ParquetStatsFileReader f(path, ParquetStatsFileReader::RecordType::NCHGCompute);
        auto first = f.begin<NCHGResult>();
        auto last = f.end<NCHGResult>();
        if (first != last) [[likely]] {
          heads.emplace_back(std::move(first));
          tails.emplace_back(std::move(last));
        }
      }
    }
  }
  return result;
}

[[nodiscard]] static FileIteratorPairs init_file_iterator_pairs(const hictk::Reference &chroms,
                                                                const MergeConfig &c) {
  try {
    const auto t0 = std::chrono::steady_clock::now();

    const auto report_path = generate_report_name(c.input_prefix);
    auto its =
        c.ignore_report_file
            ? init_file_iterator_pairs(c.input_prefix, chroms)
            : init_file_iterator_pairs(NCHGResultMetadata::from_file(report_path), c.input_prefix);

    if (its.heads.empty()) {
      throw std::runtime_error("unable to find any non-empty table");
    }

    const auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO("enumerated {} non-empty table(s) in {}", its.heads.size(),
                format_duration(t1 - t0));
    return its;
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to enumerate tables under prefix \"{}\": {}",
                                         c.input_prefix.string(), e.what()));
  }
}

using RecordQueue = moodycamel::BlockingConcurrentQueue<NCHGResult>;

[[nodiscard]] static std::size_t producer_fx(const FileIteratorPairs &iterator_pairs,
                                             RecordQueue &queue, std::atomic<bool> &early_return) {
  try {
    std::size_t records_enqueued{};
    const KMerger merger(iterator_pairs.heads, iterator_pairs.tails);

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
    ParquetStatsFileWriter writer(chromosomes, c.output_path, c.force, c.compression_method,
                                  c.compression_lvl, c.threads - 2);

    std::size_t records_dequeued = 0;
    NCHGResult buffer{};

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

      writer.append(buffer);
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

    writer.finalize<NCHGResult>();

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

static void validate_input_files(const std::filesystem::path &input_prefix) {
  const auto path_to_report = generate_report_name(input_prefix);
  SPDLOG_INFO("using \"{}\" to validate input files...", path_to_report);
  if (!std::filesystem::exists(path_to_report)) {
    throw std::runtime_error(
        fmt::format("unable to verify input file(s) integrity: file \"{}\" is missing: no such "
                    "file or directory",
                    path_to_report));
  }

  try {
    const auto t0 = std::chrono::steady_clock::now();
    const auto metadata = NCHGResultMetadata::from_file(path_to_report);
    metadata.validate();
    const auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO("SUCCESS! Validated {} files in {}", metadata.records().size(),
                format_duration(t1 - t0));
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("input files validation failed: {}", e.what()));
  }
}

[[nodiscard]] static hictk::Reference import_chromosomes(
    const std::filesystem::path &input_prefix) {
  try {
    const auto path_to_chrom_sizes = generate_chrom_sizes_name(input_prefix);
    SPDLOG_INFO("reading chromosomes from \"{}\"...", path_to_chrom_sizes);
    auto chroms = hictk::Reference::from_chrom_sizes(path_to_chrom_sizes);
    SPDLOG_INFO("read {} chromosomes!", chroms.size());
    return chroms;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to read chromosomes from \"{}\": {}", input_prefix, e.what()));
  }
}

int run_command(const MergeConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();

  if (!c.ignore_report_file) {
    validate_input_files(c.input_prefix);
  }
  const auto chroms = import_chromosomes(c.input_prefix);
  const auto iterator_pairs = init_file_iterator_pairs(chroms, c);

  moodycamel::BlockingConcurrentQueue<NCHGResult> queue(64 * 1024);
  std::atomic<bool> early_return{false};

  auto producer = std::async(std::launch::deferred, [&] {
    SPDLOG_DEBUG("spawning producer thread...");
    return producer_fx(iterator_pairs, queue, early_return);
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
