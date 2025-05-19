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

#include <fmt/format.h>
#include <moodycamel/blockingconcurrentqueue.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <future>
#include <hictk/chromosome.hpp>
#include <hictk/reference.hpp>
#include <stdexcept>
#include <string>
#include <tuple>
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

class ParquetFileMerger {
  using RecordIterator = ParquetStatsFileReader::iterator<NCHGResult>;
  using Record = std::tuple<hictk::Chromosome, hictk::Chromosome, std::filesystem::path>;

  std::queue<Record> _files;
  hictk::Chromosome _chrom1{};

 public:
  ParquetFileMerger(const hictk::Reference &chroms, std::vector<std::filesystem::path> files)
      : _files(init_file_queue(chroms, files)) {}

  auto next_chunk() -> std::optional<KMerger<RecordIterator>> {
    if (_files.empty()) {
      _chrom1 = {};
      return {};
    }

    _chrom1 = std::get<0>(_files.front());
    std::vector<std::string_view> chroms{};
    SPDLOG_DEBUG("ParquetFileMerger: processing chunk for {}", _chrom1.name());

    std::vector<RecordIterator> heads{};
    std::vector<RecordIterator> tails{};
    while (std::get<0>(_files.front()) == _chrom1) {
      const auto &[_, chrom2, path] = _files.front();
      chroms.emplace_back(chrom2.name());
      auto [first, last] = get_iterators_from_file(path);
      _files.pop();
      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
    }

    SPDLOG_DEBUG("ParquetFileMerger: merging {} chunks for {}...", heads.size(), _chrom1.name());
    return KMerger{heads, tails};
  }

 private:
  [[nodiscard]] static std::size_t find_dot_or_throw(std::string_view str, std::size_t offset) {
    const auto pos = str.rfind('.', offset);
    if (pos != std::string::npos) {
      return pos;
    }

    throw std::runtime_error(fmt::format(
        "invalid file name \"{}\": file names should be like prefix.chrA.chrB.parquet", str));
  }

  [[nodiscard]] static hictk::Chromosome fetch_chrom_or_throw(const hictk::Reference &chroms,
                                                              std::string_view chrom_name) {
    try {
      return chroms.at(chrom_name);
    } catch (const std::out_of_range &) {
      throw std::out_of_range(
          fmt::format("unable to find chromosome \"{}\": chromosome not found in "
                      "reference assembly",
                      chrom_name));
    }
  }

  static auto init_file_queue(const hictk::Reference &chroms,
                              std::vector<std::filesystem::path> &files) -> std::queue<Record> {
    std::vector<Record> records(files.size());
    std::ranges::transform(
        std::ranges::views::as_rvalue(files), records.begin(), [&](std::filesystem::path &&file) {
          const auto name = file.stem().string();

          const auto offset2 = find_dot_or_throw(name, name.size());
          const auto offset1 = find_dot_or_throw(name, offset2 - 1);

          return Record{
              fetch_chrom_or_throw(chroms, name.substr(offset1 + 1, offset2 - (offset1 + 1))),
              fetch_chrom_or_throw(chroms, name.substr(offset2 + 1)), std::move(file)};
        });

    std::ranges::sort(records);
    return {std::make_move_iterator(records.begin()), std::make_move_iterator(records.end())};
  }

  static auto get_iterators_from_file(const std::filesystem::path &path)
      -> std::pair<RecordIterator, RecordIterator> {
    try {
      ParquetStatsFileReader f(path, ParquetStatsFileReader::RecordType::NCHGCompute);
      auto first = f.begin<NCHGResult>();
      auto last = f.end<NCHGResult>();
      assert(first != last);
      return std::make_pair(std::move(first), std::move(last));
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format("failed to initialize iterators for file \"{}\": {}", path, e.what()));
    }
  }
};

[[nodiscard]] static std::string generate_chrom_sizes_name(
    const std::filesystem::path &input_prefix) {
  return fmt::format("{}.chrom.sizes", input_prefix.string());
}

[[nodiscard]] static std::string generate_report_name(const std::filesystem::path &input_prefix) {
  return fmt::format("{}.json", input_prefix);
}

[[nodiscard]] static std::vector<std::filesystem::path> enumerate_parquet_tables(
    const NCHGResultMetadata &metadata, const std::filesystem::path &input_prefix) {
  assert(!metadata.records().empty());

  SPDLOG_INFO("enumerating tables based on file(s) listed in report file \"{}\"...",
              generate_report_name(input_prefix));

  const auto root_folder =
      input_prefix.has_parent_path() ? input_prefix.parent_path() : std::filesystem::current_path();

  std::vector<std::filesystem::path> files;
  files.reserve(metadata.records().size() - 1);

  for (const auto &record : metadata.records()) {
    if (record.name.extension() != ".parquet") {
      continue;
    }

    auto path = root_folder / record.name;
    ParquetStatsFileReader f(path, ParquetStatsFileReader::RecordType::NCHGCompute);
    if (f.begin<NCHGResult>() != f.end<NCHGResult>()) {
      files.emplace_back(std::move(path));
    }
  }

  return files;
}

[[nodiscard]] static std::vector<std::filesystem::path> enumerate_parquet_tables(
    const std::filesystem::path &input_prefix, const hictk::Reference &chroms) {
  SPDLOG_INFO("enumerating table based on chromosome(s) listed in \"{}\"...",
              generate_chrom_sizes_name(input_prefix));
  std::vector<std::filesystem::path> files;

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

      std::filesystem::path path =
          fmt::format("{}.{}.{}.parquet", input_prefix.string(), chrom1.name(), chrom2.name());
      if (std::filesystem::exists(path)) {
        ParquetStatsFileReader f(path, ParquetStatsFileReader::RecordType::NCHGCompute);
        if (f.begin<NCHGResult>() != f.end<NCHGResult>()) {
          files.emplace_back(std::move(path));
        }
      }
    }
  }
  return files;
}

[[nodiscard]] static std::vector<std::filesystem::path> enumerate_parquet_tables(
    const hictk::Reference &chroms, const MergeConfig &c) {
  try {
    const auto t0 = std::chrono::steady_clock::now();

    const auto report_path = generate_report_name(c.input_prefix);
    auto files =
        c.ignore_report_file
            ? enumerate_parquet_tables(c.input_prefix, chroms)
            : enumerate_parquet_tables(NCHGResultMetadata::from_file(report_path), c.input_prefix);

    if (files.empty()) {
      throw std::runtime_error("unable to find any non-empty table");
    }

    const auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO("enumerated {} non-empty table(s) in {}", files.size(), format_duration(t1 - t0));
    return files;
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format("failed to enumerate tables under prefix \"{}\": {}",
                                         c.input_prefix.string(), e.what()));
  }
}

using RecordQueue = moodycamel::BlockingConcurrentQueue<NCHGResult>;

[[nodiscard]] static std::size_t producer_fx(
    const hictk::Reference &chroms, const std::vector<std::filesystem::path> &parquet_files,
    RecordQueue &queue, std::atomic<bool> &early_return) {
  try {
    std::size_t records_enqueued{};
    ParquetFileMerger file_merger(chroms, parquet_files);

    while (true) {
      auto chunk_merger = file_merger.next_chunk();
      if (!chunk_merger.has_value()) {
        break;
      }

      for (const auto &s : *chunk_merger) {
        while (!queue.try_enqueue(s)) [[unlikely]] {
          if (early_return) [[unlikely]] {
            return records_enqueued;
          }
        }
        ++records_enqueued;
      }
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

[[nodiscard]] static constexpr spdlog::level::level_enum extract_log_lvl(std::string_view msg) {
  if (msg.contains(" [critical]: ")) [[unlikely]] {
    return spdlog::level::critical;
  }
  if (msg.contains(" [error]: ")) [[unlikely]] {
    return spdlog::level::err;
  }
  if (msg.contains(" [warning]: ")) [[unlikely]] {
    return spdlog::level::warn;
  }
  return spdlog::level::info;
}

[[nodiscard]] static constexpr bool is_warning_replay_header(std::string_view msg) {
  auto offset = msg.find(" [warning]: ");
  if (offset == std::string_view::npos) [[likely]] {
    return false;
  }
  offset = msg.find("replaying the last ", offset);
  if (offset == std::string_view::npos) {
    return false;
  }
  offset = msg.find("warning message(s)", offset);
  return offset != std::string_view::npos;
}

[[nodiscard]] static std::vector<std::string> try_parse_nchg_compute_log_file(
    const std::filesystem::path &path) {
  if (!std::filesystem::exists(path)) {
    return {};
  }

  std::vector<std::string> messages{};
  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  try {
    fs.open(path);
    std::string buffer{};
    auto min_log_level = spdlog::level::warn;

    for (std::size_t i = 1; std::getline(fs, buffer); ++i) {
      const auto log_level = extract_log_lvl(buffer);
      if (log_level >= min_log_level && is_warning_replay_header(buffer)) {
        min_log_level = spdlog::level::err;
      }

      if (log_level >= min_log_level) {
        messages.emplace_back(fmt::format("line {}: {}", i, buffer));
      }
    }
  } catch (...) {
    if (fs.eof()) {
      return messages;
    }
  }
  return {};
}

[[noreturn]] static void handle_incomplete_output(const std::filesystem::path &report_file) {
  constexpr auto msg =
      "input files validation failed: detected incomplete and/or corrupted input file(s).\n"
      "This suggests that a previous run of NCHG compute did not complete successfully.\n"
      "Possible reasons are:\n"
      " - NCHG compute encountered an error\n"
      " - NCHG compute was killed by the OS due to e.g. running out of memory\n"
      " - NCHG compute was killed by a workflow manager, such as SLURM, due to e.g. running "
      "out of time or resources\n"
      "\n"
      "If you are sure that all files produced by NCHG compute are available and are not "
      "corrupted, you can pass --ignore-report-file to NCHG merge to disable input file "
      "validation.";

  auto log_file_name = report_file;
  log_file_name.replace_extension(".log");
  if (!std::filesystem::exists(log_file_name)) {
    throw std::runtime_error(msg);
  }
  const auto log_entries_of_interest = try_parse_nchg_compute_log_file(log_file_name);

  if (log_entries_of_interest.empty()) {
    throw std::runtime_error(
        fmt::format("{}\n"
                    "An attempt was made to parse log file \"{}\" but the file does not seem to "
                    "contain helpful information.",
                    msg, log_file_name.string()));
  }

  throw std::runtime_error(
      fmt::format("{}\n"
                  "Parsing log file \"{}\" identified the following log messages as potentially "
                  "useful for further troubleshooting:\n{}",
                  msg, log_file_name.string(), fmt::join(log_entries_of_interest, "\n")));
}

[[noreturn]] static void handle_invalid_report_file(
    const std::filesystem::path &path, const NCHGResultMetadata::ValidationResult &res) {
  try {
    res.throw_exception();
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(fmt::format(
        "input files validation failed: report file \"{}\" is corrupted or incomplete:\n{}",
        path.string(), e.what()));
  }
}

[[nodiscard]] static auto classify_invalid_records(
    const NCHGResultMetadata::ValidationResult &res) {
  struct InvalidRecordsReport {
    std::vector<std::pair<std::filesystem::path, std::string>> missing_files{};
    std::vector<std::pair<std::filesystem::path, std::string>> size_mismatch{};
    std::vector<std::pair<std::filesystem::path, std::string>> checksum_mismatch{};
    std::vector<std::pair<std::filesystem::path, std::string>> other_errors{};
  };

  InvalidRecordsReport result{};
  for (const auto &[file, msg] : res.record_validation_failures) {
    if (msg.contains("file does not exist")) {
      result.missing_files.emplace_back(file, msg);
    } else if (msg.contains("file size mismatch: ")) {
      constexpr std::string_view query{"file size mismatch: "};
      const auto offset = msg.find(query) + query.size();
      result.size_mismatch.emplace_back(file, msg.substr(offset));
    } else if (msg.contains("checksum mismatch: ")) {
      constexpr std::string_view query{"checksum mismatch: "};
      const auto offset = msg.find(query) + query.size();
      result.checksum_mismatch.emplace_back(file, msg.substr(offset));
    } else {
      result.other_errors.emplace_back(file, msg);
    }
  }

  return result;
}

[[nodiscard]] static std::string generate_suggestions(
    const std::filesystem::path &path,
    const std::vector<std::pair<std::filesystem::path, std::string>> &missing_files,
    const std::vector<std::pair<std::filesystem::path, std::string>> &size_mismatch,
    const std::vector<std::pair<std::filesystem::path, std::string>> &checksum_mismatch) {
  if (missing_files.empty() && size_mismatch.empty() && checksum_mismatch.empty()) {
    return "";
  }

  std::string suggestions{"Suggestions:"};
  if (!missing_files.empty()) {
    suggestions += fmt::format(
        "\n - Make sure no file has been deleted or renamed.\n"
        "   If files have been move to a different folder, please make sure that all files "
        "listed in \"{}\" have been moved to the new location.",
        std::filesystem::canonical(path).string());
  }

  if (size_mismatch.size() + checksum_mismatch.size() != 0) {
    suggestions +=
        "\n - Make sure that files with a different size or checksum have not been overwritten "
        "or modified after their creation.\n"
        "   If you are sure that all files produced by NCHG compute are available and are not "
        "corrupted, you can pass --ignore-report-file to NCHG merge to disable input file "
        "validation.";
  }

  return suggestions;
}

[[nodiscard]] static std::string report_missing_files(
    const std::vector<std::pair<std::filesystem::path, std::string>> &missing_files) {
  assert(!missing_files.empty());
  constexpr std::string_view prefix = "\n    - ";
  std::vector<std::string> buffer{};
  std::ranges::transform(missing_files, std::back_inserter(buffer),
                         [](const auto &kv) { return kv.first.string(); });

  return fmt::format("\n - Missing file(s):{}{}", prefix, fmt::join(buffer, prefix));
}

[[nodiscard]] static std::string report_size_mismatch(
    const std::vector<std::pair<std::filesystem::path, std::string>> &size_mismatch) {
  assert(!size_mismatch.empty());
  constexpr std::string_view prefix = "\n    - ";
  std::vector<std::string> buffer{};
  std::ranges::transform(size_mismatch, std::back_inserter(buffer), [](const auto &kv) {
    return fmt::format("{}: {}", kv.first.string(), kv.second);
  });
  return fmt::format("\n - File size mismatch:{}{}", prefix, fmt::join(buffer, prefix));
}

[[nodiscard]] static std::string report_checksum_mismatch(
    const std::vector<std::pair<std::filesystem::path, std::string>> &checksum_mismatch) {
  assert(!checksum_mismatch.empty());
  constexpr std::string_view prefix = "\n    - ";
  std::vector<std::string> buffer{};
  std::ranges::transform(checksum_mismatch, std::back_inserter(buffer), [](const auto &kv) {
    return fmt::format("{}: {}", kv.first.string(), kv.second);
  });

  return fmt::format("\n - File checksum mismatch:{}{}", prefix, fmt::join(buffer, prefix));
}

[[nodiscard]] static std::string report_other_errors(
    const std::vector<std::pair<std::filesystem::path, std::string>> &other_errors) {
  assert(!other_errors.empty());
  constexpr std::string_view prefix = "\n    - ";
  std::vector<std::string> buffer{};
  std::ranges::transform(other_errors, std::back_inserter(buffer), [](const auto &kv) {
    return fmt::format("{}: {}", kv.first.string(), kv.second);
  });

  return fmt::format("\n - Other errors:{}{}", prefix, fmt::join(buffer, prefix));
}

[[noreturn]] static void handle_invalid_records(const std::filesystem::path &path,
                                                const NCHGResultMetadata::ValidationResult &res) {
  const auto [missing_files, size_mismatch, checksum_mismatch, other_errors] =
      classify_invalid_records(res);

  std::string msg{"input files validation failed:"};

  if (!missing_files.empty()) {
    msg += report_missing_files(missing_files);
  }

  if (!size_mismatch.empty()) {
    msg += report_size_mismatch(size_mismatch);
  }

  if (!checksum_mismatch.empty()) {
    msg += report_checksum_mismatch(checksum_mismatch);
  }

  if (!other_errors.empty()) {
    msg += report_other_errors(other_errors);
  }

  const auto suggestions =
      generate_suggestions(path, missing_files, size_mismatch, checksum_mismatch);

  std::string summary{"Summary:"};
  if (!missing_files.empty()) {
    summary += fmt::format("\n - {} record(s) refer to non-existing files", missing_files.size());
  }
  if (!size_mismatch.empty() || !checksum_mismatch.empty()) {
    summary += fmt::format("\n - {} file(s) have changed after their creation",
                           size_mismatch.size() + checksum_mismatch.size());
  }

  throw std::runtime_error(
      fmt::format("{}{}{}\n{}", msg, suggestions.empty() ? "" : "\n", suggestions, summary));
}

static void handle_validation_errors(const std::filesystem::path &path,
                                     const NCHGResultMetadata::ValidationResult &res) {
  if (!!res) [[likely]] {
    return;
  }

  if (!res.successfully_finalized) {
    handle_incomplete_output(path);
  }

  if (!res.report_validation_failures.empty()) {
    handle_invalid_report_file(path, res);
  }

  if (!res.record_validation_failures.empty()) {
    handle_invalid_records(path, res);
  }

  try {
    res.throw_exception();
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(
        fmt::format("input files validation failed: an unhandled exception occurred while parsing "
                    "report file \"{}\": {}",
                    path.string(), e.what()));
  }
}

static void validate_input_files(const std::filesystem::path &input_prefix, std::size_t threads) {
  const auto path_to_report = generate_report_name(input_prefix);
  SPDLOG_INFO("using \"{}\" to validate input files...", path_to_report);
  if (!std::filesystem::exists(path_to_report)) {
    throw std::runtime_error(
        fmt::format("unable to verify input file(s) integrity: file \"{}\" is missing: no such "
                    "file or directory",
                    path_to_report));
  }

  const auto t0 = std::chrono::steady_clock::now();
  const auto metadata = NCHGResultMetadata::from_file(path_to_report, false);
  std::unique_ptr<BS::light_thread_pool> tpool =
      threads > 1 ? std::make_unique<BS::light_thread_pool>(threads) : nullptr;
  const auto validation_result = metadata.validate(tpool.get());
  handle_validation_errors(path_to_report, validation_result);
  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("SUCCESS! Validated {} files in {}", metadata.records().size(),
              format_duration(t1 - t0));
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
    validate_input_files(c.input_prefix, c.threads);
  }
  const auto chroms = import_chromosomes(c.input_prefix);
  const auto parquet_files = enumerate_parquet_tables(chroms, c);

  moodycamel::BlockingConcurrentQueue<NCHGResult> queue(64 * 1024);
  std::atomic<bool> early_return{false};

  auto producer = std::async(std::launch::deferred, [&] {
    SPDLOG_DEBUG("spawning producer thread...");
    return producer_fx(chroms, parquet_files, queue, early_return);
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
