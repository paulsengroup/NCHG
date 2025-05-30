// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
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
// <https://www.gnu.org/licenses/>

#include "./scheduler.hpp"

#include <arrow/io/file.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <parallel_hashmap/btree.h>
#include <parquet/stream_reader.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <atomic>
#include <boost/asio/io_context.hpp>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

#include "./genomic_domains.hpp"
#include "./logging.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/file_store.hpp"
#include "nchg/genomic_domains.hpp"
#include "nchg/nchg.hpp"
#include "nchg/parquet_stats_file_writer.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tmpdir.hpp"

// clang-format off
// As of HighFive 2.9.0, these headers must be included after HighFive/hictk,
// otherwise this source file fails to compile with MSVC
#include <boost/process/v2/environment.hpp>
#include <boost/process/v2/process.hpp>
#include <boost/process/v2/stdio.hpp>
// clang-format on

namespace nchg {
class ProcessContext {
  using LockedContext = std::pair<std::unique_lock<std::mutex>, boost::asio::io_context *>;
  boost::asio::io_context _ctx;
  std::mutex _mtx;

 public:
  ProcessContext() = default;
  [[nodiscard]] LockedContext operator()() { return {std::unique_lock{_mtx}, &_ctx}; }
};

[[nodiscard]] static auto init_nchg(const std::shared_ptr<const hictk::File> &f,
                                    const std::optional<ExpectedValues> &expected_values,
                                    const ComputePvalConfig &c) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());

  struct Result {
    NCHG nchg;
    std::shared_ptr<const std::vector<bool>> bin1_mask;
    std::shared_ptr<const std::vector<bool>> bin2_mask;
  };

  const auto &chrom1 = f->chromosomes().at(*c.chrom1);
  const auto &chrom2 = f->chromosomes().at(*c.chrom2);

  if (expected_values.has_value()) {
    auto [bin1_mask, bin2_mask] = expected_values->bin_mask(chrom1, chrom2);
    return Result{.nchg = NCHG{f, chrom1, chrom2, *expected_values},
                  .bin1_mask = std::move(bin1_mask),
                  .bin2_mask = std::move(bin2_mask)};
  }

  if (!c.path_to_expected_values.empty()) {
    SPDLOG_INFO("reading expected values from {}...", c.path_to_expected_values);
    return init_nchg(f, ExpectedValues::deserialize(c.path_to_expected_values), c);
  }

  const auto bin_mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  const auto evs = ExpectedValues::chromosome_pair(
      f, chrom1, chrom2,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      bin_mask);

  return init_nchg(f, evs, c);
}

[[nodiscard]] static std::size_t process_domains(
    const std::shared_ptr<const hictk::File> &f, const GenomicDomains &domains,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  assert(f);
  assert(std::filesystem::exists(c.path_to_domains));
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(!c.output_path.empty());

  SPDLOG_INFO("[{}:{}]: begin processing domains from {}...", *c.chrom1, *c.chrom2,
              c.path_to_domains);

  ParquetStatsFileWriter writer(f->chromosomes(), c.output_path, c.force, c.compression_method,
                                c.compression_lvl, c.threads);

  const auto &chrom1 = f->chromosomes().at(*c.chrom1);
  const auto &chrom2 = f->chromosomes().at(*c.chrom2);

  if (!domains.contains(chrom1, chrom2)) {
    writer.finalize<NCHGResult>();
    return 0;
  }

  const auto [nchg, bin1_mask, bin2_mask] = init_nchg(f, expected_values, c);
  assert(bin1_mask);
  assert(bin2_mask);

  const auto domains_with_interactions = map_interactions_to_domains(
      *f, domains, nchg.expected_matrix(), chrom1, chrom2, c.min_delta, c.max_delta, *bin1_mask,
      *bin2_mask, c.domain_aggregation_stategy);

  std::size_t num_records = 0;
  for (const auto &[domain, obs, exp] : domains_with_interactions) {
    const auto s = nchg.compute(domain, obs, exp, c.bad_bin_fraction);

    if (std::isfinite(s.odds_ratio) && s.omega != 0) [[likely]] {
      writer.append(s);
    }
    ++num_records;
  }

  writer.finalize<NCHGResult>();
  return num_records;
}

[[nodiscard]] static std::size_t process_one_chromosome_pair(
    const std::shared_ptr<const hictk::File> &f,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(!c.output_path.empty());

  SPDLOG_INFO("[{}:{}]: begin processing interactions...", *c.chrom1, *c.chrom2);

  const auto &chrom1 = f->chromosomes().at(*c.chrom1);
  const auto &chrom2 = f->chromosomes().at(*c.chrom2);
  const auto nchg = init_nchg(f, expected_values, c).nchg;

  ParquetStatsFileWriter writer(f->chromosomes(), c.output_path, c.force, c.compression_method,
                                c.compression_lvl, c.threads);

  std::size_t num_records = 0;
  auto first = nchg.begin(chrom1, chrom2);
  auto last = nchg.end(chrom1, chrom2);

  std::visit(
      [&]<typename It>(It &it1) {
        auto &it2 = std::get<It>(last);
        std::for_each(it1, it2, [&](const auto &s) {
          ++num_records;

          if (std::isfinite(s.odds_ratio) && s.omega != 0) [[likely]] {
            writer.append(s);
          }
        });
      },
      first);

  writer.finalize<NCHGResult>();
  return num_records;
}

[[nodiscard]] static std::string aggregation_strategy_to_str(DomainAggregationStrategy strategy) {
  switch (strategy) {
    case DomainAggregationStrategy::AUTO:
      return "auto";
    case DomainAggregationStrategy::SINGLE_PASS:
      return "one-pass";
    case DomainAggregationStrategy::MULTI_PASS:
      return "multi-pass";
  }
  unreachable_code();
}

[[nodiscard]] static boost::process::process_environment generate_subprocess_env(
    const MessageQueue &msg_queue) {
  constexpr boost::string_view queue_name_env_variable{"NCHG_LOG_MESSAGE_QUEUE_NAME"};

  static std::once_flag flag;
  using EnvironmentKV =
      std::pair<boost::process::environment::key, boost::process::environment::value>;
  static std::vector<EnvironmentKV> vars{};

  std::call_once(flag, [&] {
    if constexpr (ndebug_not_defined()) {
      std::ranges::transform(
          boost::process::environment::current(), std::back_inserter(vars), [](const auto &kv) {
            return std::make_pair(boost::process::environment::key{kv.key()},
                                  boost::process::environment::value{kv.value()});
          });
    } else {
      // NOLINTNEXTLINE(*-mt-unsafe)
      if (const auto *var = std::getenv("NCHG_CI"); var) {
        vars.emplace_back("NCHG_CI", var);
      }
      // NOLINTNEXTLINE(*-mt-unsafe)
      if (const auto *var = std::getenv("LLVM_PROFILE_FILE"); var) {
        vars.emplace_back("LLVM_PROFILE_FILE", var);
      }
    }
    vars.emplace_back(boost::process::environment::key{queue_name_env_variable},
                      boost::process::environment::value{msg_queue.name()});
  });

  return {vars};
}

[[nodiscard]] static std::vector<std::string> generate_compute_args(const hictk::Chromosome &chrom1,
                                                                    const hictk::Chromosome &chrom2,
                                                                    const ComputePvalConfig &c) {
  assert(c.output_prefix.empty());
  assert(!c.output_path.empty());

  std::vector<std::string> args{
      "compute",
      c.path_to_hic.string(),
      c.output_path.string(),
      "--chrom1",
      std::string{chrom1.name()},
      "--chrom2",
      std::string{chrom2.name()},
      "--threads",
      "1",
      "--bad-bin-fraction",
      fmt::to_string(c.bad_bin_fraction),
      "--compression-level",
      fmt::to_string(c.compression_lvl),
      "--compression-method",
      c.compression_method,
      c.skip_empty_matrices ? "--skip-empty-matrices" : "--keep-empty-matrices",
  };

  if (c.resolution.has_value()) {
    args.emplace_back("--resolution");
    args.emplace_back(fmt::to_string(*c.resolution));
  }

  if (!c.path_to_domains.empty()) {
    args.emplace_back("--domains");
    args.emplace_back(c.path_to_domains.string());
    args.emplace_back("--interaction-aggregation-strategy");
    args.emplace_back(aggregation_strategy_to_str(c.domain_aggregation_stategy));
  }

  if (c.path_to_expected_values.empty()) {
    args.emplace_back("--min-delta");
    args.emplace_back(fmt::to_string(c.min_delta));
    args.emplace_back("--max-delta");
    args.emplace_back(fmt::to_string(c.max_delta));
    args.emplace_back("--bin-aggregation-possible-distances-cutoff");
    args.emplace_back(fmt::to_string(c.bin_aggregation_possible_distances_cutoff));
    args.emplace_back("--bin-aggregation-observed-distances-cutoff");
    args.emplace_back(fmt::to_string(c.bin_aggregation_observed_distances_cutoff));
    if (c.interpolate_expected_values) {
      args.emplace_back("--interpolate-expected-values");
    } else {
      args.emplace_back("--no-interpolate-expected-values");
    }
    args.emplace_back("--evs-interpolation-qtile");
    args.emplace_back(fmt::to_string(c.interpolation_qtile));
    args.emplace_back("--evs-interpolation-window");
    args.emplace_back(fmt::to_string(c.interpolation_window_size));
    args.emplace_back("--mad-max");
    args.emplace_back(fmt::to_string(c.mad_max));
    if (!c.path_to_bin_mask.empty()) {
      args.emplace_back("--bin-mask");
      args.emplace_back(c.path_to_bin_mask.string());
    }
  } else {
    args.emplace_back("--expected-values");
    args.emplace_back(c.path_to_expected_values.string());
  }

  if (c.force) {
    args.emplace_back("--force");
  }

  return args;
}

[[nodiscard]] static std::chrono::milliseconds generate_random_sleep_time_ms(
    double target_sleep_ms, double stddev = 1.0, double min = 100.0, double max = 10'000.0) {
  assert(target_sleep_ms > 0);
  assert(max >= min);
  assert(stddev >= 0);
  std::random_device rd{};
  std::mt19937 rand_eng{rd()};
  const auto sleep_time_ms =
      std::clamp(min, max, std::normal_distribution{target_sleep_ms, stddev}(rand_eng));
  return std::chrono::milliseconds{static_cast<int>(sleep_time_ms)};
}

[[nodiscard]] static boost::process::process spawn_compute_process(
    ProcessContext &ctx, const MessageQueue &msg_queue, const ComputePvalConfig &c,
    const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) {
  SPDLOG_INFO("[{}:{}]: begin processing...", chrom1.name(), chrom2.name());

  const auto args = generate_compute_args(chrom1, chrom2, c);
  auto env = generate_subprocess_env(msg_queue);

  phmap::btree_set<std::string> errors{};
  double mean_sleep_time = 250;

  for (std::size_t attempt = 0; attempt < 10; ++attempt) {
    try {
      const auto [lck, asio_ctx] = ctx();
      boost::process::process proc(
          *asio_ctx, c.exec.string(), args, env,
          boost::process::process_stdio{.in = nullptr, .out = nullptr, .err = {}});
      if (proc.running() || proc.exit_code() == 0) {
        SPDLOG_DEBUG("[{}:{}]: spawned worker process {}...", chrom1.name(), chrom2.name(),
                     proc.id());
        return proc;
      }
      SPDLOG_WARN("[{}:{}]: spawning worker process {} failed (attempt {}/10)...", chrom1.name(),
                  chrom2.name(), proc.id(), attempt + 1);
      proc.terminate();
    } catch (const std::exception &e) {
      SPDLOG_WARN("[{}:{}]: spawning worker process failed (attempt {}/10): {}", chrom1.name(),
                  chrom2.name(), attempt + 1, e.what());
      errors.emplace(e.what());
    }

    const auto sleep_time = generate_random_sleep_time_ms(mean_sleep_time);
    SPDLOG_DEBUG("[{}:{}]: sleeping {} before attempting to spawn worker process one more time...",
                 chrom1.name(), chrom2.name(), sleep_time);
    std::this_thread::sleep_for(sleep_time);
    mean_sleep_time *= 1.5;
  }

  if (errors.empty()) {
    throw std::runtime_error(fmt::format("failed to spawn worker process: {} {}", c.exec.string(),
                                         fmt::join(args, " ")));
  }
  if (errors.size() == 1) {
    throw std::runtime_error(fmt::format("failed to spawn worker process: {} {}: {}",
                                         c.exec.string(), fmt::join(args, " "), *errors.begin()));
  }
  throw std::runtime_error(
      fmt::format("failed to spawn worker process: {} {}\n"
                  "Exception(s):\n - {}",
                  c.exec.string(), fmt::join(args, " "), fmt::join(errors, "\n - ")));
}

[[nodiscard]] static std::size_t worker_fx(FileStore &file_store, const hictk::Chromosome &chrom1,
                                           const hictk::Chromosome &chrom2,
                                           const std::optional<GenomicDomains> &domains,
                                           const TmpDir &tmpdir, const MessageQueue &msg_queue,
                                           ProcessContext &ctx, const ComputePvalConfig &config,
                                           bool trans_expected_values_avail,
                                           std::atomic<bool> &early_return) {
  if (early_return) {
    return 0;
  }

  try {
    const auto t0 = std::chrono::steady_clock::now();
    auto child_config = config;

    child_config.output_path =
        generate_output_file_name(config.output_prefix.string(), chrom1, chrom2);
    child_config.output_prefix.clear();

    std::filesystem::path domain_file{};
    if (domains.has_value()) {
      domain_file = write_domains_to_file(*domains, tmpdir(), chrom1, chrom2, child_config.force);
      child_config.path_to_domains = domain_file;
    } else {
      assert(child_config.path_to_domains.empty());
    }

    if (!trans_expected_values_avail && chrom1 != chrom2) {
      child_config.path_to_expected_values.clear();
    }

    auto proc = spawn_compute_process(ctx, msg_queue, child_config, chrom1, chrom2);
    try {
      proc.wait();
    } catch (const std::exception &e) {
      const std::string_view msg{e.what()};
      // Deal with processes that terminated almost instantaneously
      if (msg.find("wait failed: No child processes") == std::string_view::npos) {
        throw;
      }
    }

    if (!domain_file.empty()) {
      std::filesystem::remove(domain_file);  // NOLINT
    }

    if (proc.exit_code() != 0) {
      early_return = true;
      throw std::runtime_error(
          fmt::format("child process terminated with code {}", proc.exit_code()));
    }
    SPDLOG_DEBUG("[{}:{}]: worker process returned with exit code 0", chrom1.name(), chrom2.name());

    file_store.register_file(child_config.output_path);

    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(child_config.output_path));
    const auto records_processed =
        parquet::StreamReader{parquet::ParquetFileReader::Open(fp)}.num_rows();

    const auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO("[{}:{}]: processed {} records in {}", chrom1.name(), chrom2.name(),
                records_processed, format_duration(t1 - t0));
    return static_cast<std::size_t>(records_processed);

  } catch (const std::exception &e) {
    early_return = true;
    throw std::runtime_error(fmt::format("error in the worker thread processing {}:{}: {}",
                                         chrom1.name(), chrom2.name(), e.what()));
  } catch (...) {
    SPDLOG_ERROR("an unknown error occurred in worker thread processing {}:{}", chrom1.name(),
                 chrom2.name());
    early_return = true;
    throw;
  }
}

[[nodiscard]] static ComputePvalConfig init_base_config(
    const ComputePvalConfig &c, const TmpDir &tmpdir,
    const std::optional<ExpectedValues> &expected_values) {
  auto base_config = c;
  if (c.path_to_expected_values.empty() && expected_values.has_value()) {
    assert(!base_config.output_prefix.empty());
    base_config.path_to_expected_values =
        tmpdir() / fmt::format("{}_expected_values.h5", base_config.output_prefix.stem().string());

    if (!base_config.force && std::filesystem::exists(base_config.path_to_expected_values)) {
      throw std::runtime_error(
          fmt::format("Refusing to overwrite file {}. Pass --force to overwrite.",
                      base_config.path_to_expected_values));
    }

    std::filesystem::remove(base_config.path_to_expected_values);  // NOLINT
    expected_values->serialize(base_config.path_to_expected_values);
  }

  return base_config;
}

std::size_t process_chromosome_pair(const ComputePvalConfig &c,
                                    const std::optional<GenomicDomains> &domains,
                                    const std::optional<ExpectedValues> &expected_values) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());

  const auto f = std::make_shared<const hictk::File>(c.path_to_hic.string(), c.resolution);

  if (domains.has_value()) {
    return process_domains(f, *domains, expected_values, c);
  }

  if (!c.path_to_domains.empty()) {
    const auto chrom1 = std::make_optional(f->chromosomes().at(*c.chrom1));
    const auto chrom2 = std::make_optional(f->chromosomes().at(*c.chrom2));

    const auto domains_ =
        parse_domains(f->chromosomes(), c.path_to_domains, true, true, chrom1, chrom2);
    return process_domains(f, domains_, expected_values, c);
  }

  return process_one_chromosome_pair(f, expected_values, c);
}

static void process_log_messages(MessageQueue &msg_queue, std::atomic<bool> &early_return) {
  std::size_t num_except = 0;
  SPDLOG_DEBUG("starting logger thread...");

  for ([[maybe_unused]] std::size_t i = 0; !early_return; ++i) {
    try {
      if (!msg_queue.receive()) {
        SPDLOG_DEBUG("logger thread is returning: processed a total of {} records", i);
        return;
      }
    } catch (const std::exception &e) {
      if (++num_except > 10) {
        early_return = true;
        throw std::runtime_error(
            fmt::format("logger thread is encountered the following exception: {}", e.what()));
      }
      SPDLOG_WARN("logger thread is encountered the following exception: {}", e.what());
    } catch (...) {
      if (++num_except > 10) {
        early_return = true;
        throw std::runtime_error("logger thread is encountered an unknown exception");
      }
      SPDLOG_WARN("logger thread is encountered an unknown exception");
    }
  }

  if (early_return) {
    SPDLOG_DEBUG("logger thread: early return signal received: returning immediately!");
  }
}

static std::size_t process_queries_mt(BS::light_thread_pool &tpool, FileStore &file_store,
                                      const ChromosomePairs &chrom_pairs,
                                      const std::optional<GenomicDomains> &domains,
                                      const std::optional<ExpectedValues> &expected_values,
                                      const TmpDir &tmpdir, const ComputePvalConfig &c) {
  const auto user_provided_expected_values = !c.path_to_expected_values.empty();
  const auto base_config = init_base_config(c, tmpdir, expected_values);

  ProcessContext ctx;
  std::atomic early_return{false};
  BS::multi_future<std::size_t> workers;
  workers.reserve(chrom_pairs.size());

  auto msg_queue = MessageQueue::create(TmpDir::generate_random_file_name("nchg-log-queue-", 32),
                                        tpool.get_thread_count(), early_return);

  auto logger = tpool.submit_task([&] { process_log_messages(msg_queue, early_return); });

  auto it = chrom_pairs.begin();
  std::chrono::milliseconds sleep_time{2};
  for (std::size_t tasks_submitted = 0; tasks_submitted < workers.capacity() && !early_return;) {
    if (tpool.get_tasks_running() == tpool.get_thread_count()) {
      std::this_thread::sleep_for(sleep_time);
      sleep_time = std::min(std::chrono::milliseconds{500}, sleep_time * 2);
      continue;
    }
    workers.emplace_back(tpool.submit_task([&, tasks_submitted, chrom_pair = *it++] {
      std::ignore = tasks_submitted;
      const auto &[chrom1, chrom2] = chrom_pair;
      SPDLOG_DEBUG("submitting task {}/{} ({}:{})...", tasks_submitted + 1, workers.size(),
                   chrom1.name(), chrom2.name());
      // NOLINTBEGIN(clang-analyzer-unix.BlockInCriticalSection)
      return worker_fx(file_store, chrom1, chrom2, domains, tmpdir, msg_queue, ctx, base_config,
                       user_provided_expected_values, early_return);
      // NOLINTEND(clang-analyzer-unix.BlockInCriticalSection)
    }));
    ++tasks_submitted;
    sleep_time = std::chrono::milliseconds{2};
  }

  std::size_t num_records = 0;
  for (const auto result : workers.get()) {
    num_records += result;
  }
  msg_queue.send_eoq_signal();
  logger.get();

  return num_records;
}

static std::size_t process_queries_st(FileStore &file_store, const ChromosomePairs &chrom_pairs,
                                      const std::optional<GenomicDomains> &domains,
                                      const std::optional<ExpectedValues> &expected_values,
                                      const ComputePvalConfig &c) {
  assert(!c.output_prefix.empty());
  std::size_t tot_num_records = 0;
  for (const auto &[chrom1, chrom2] : chrom_pairs) {
    try {
      auto config = c;
      config.chrom1 = chrom1.name();
      config.chrom2 = chrom2.name();
      config.compute_cis = true;
      config.compute_trans = true;

      config.output_path =
          fmt::format("{}.{}.{}.parquet", c.output_prefix.string(), chrom1.name(), chrom2.name());
      config.output_prefix.clear();

      const auto t0 = std::chrono::steady_clock::now();
      if (config.chrom1 == config.chrom2) {
        assert(expected_values.has_value());
      }
      const auto num_records = config.chrom1 == config.chrom2
                                   ? process_chromosome_pair(config, domains, expected_values)
                                   : process_chromosome_pair(config, domains);
      tot_num_records += num_records;

      file_store.register_file(config.output_path);

      const auto t1 = std::chrono::steady_clock::now();
      SPDLOG_INFO("[{}:{}]: processed {} records in {}!", chrom1.name(), chrom2.name(), num_records,
                  format_duration(t1 - t0));
      SPDLOG_INFO("[{}:{}]: {} records have been written to file \"{}\"", chrom1.name(),
                  chrom2.name(), num_records, config.output_path.string());

    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format("error in while processing {}:{}: {}", chrom1.name(),
                                           chrom2.name(), e.what()));
    } catch (...) {
      throw std::runtime_error(fmt::format("An unknown error occurred while processing {}:{}",
                                           chrom1.name(), chrom2.name()));
    }
  }

  return tot_num_records;
}

std::size_t dispatch_queries(FileStore &file_store, const ChromosomePairs &chrom_pairs,
                             const std::optional<GenomicDomains> &domains,
                             const std::optional<ExpectedValues> &expected_values,
                             const ComputePvalConfig &c) {
  const std::filesystem::path chrom_sizes_path{generate_chrom_sizes_file_name(c.output_prefix)};
  write_chrom_sizes_to_file(hictk::File(c.path_to_hic.string(), c.resolution).chromosomes(),
                            chrom_sizes_path, c.force);
  file_store.register_file(chrom_sizes_path);

  if (c.threads > 1) {
    assert(!c.output_prefix.parent_path().empty());
    const TmpDir tmpdir{c.output_prefix.parent_path(), true};
    SPDLOG_INFO("writing temporary files under folder \"{}\"...", tmpdir());

    auto num_threads = c.threads;
    if (c.threads > chrom_pairs.size()) {
      num_threads = chrom_pairs.size();
      SPDLOG_WARN(
          "number of threads specified through --threads exceeds the number of chromosome pairs "
          "to be processed: limiting concurrency to {} thread(s)",
          num_threads);
    }
    assert(num_threads != 0);
    BS::light_thread_pool tpool(num_threads + 1);
    return process_queries_mt(tpool, file_store, chrom_pairs, domains, expected_values, tmpdir, c);
  }
  return process_queries_st(file_store, chrom_pairs, domains, expected_values, c);
}

std::filesystem::path generate_report_file_name(const std::filesystem::path &output_prefix) {
  assert(!output_prefix.empty());
  return fmt::format("{}.json", output_prefix);
}

std::filesystem::path generate_chrom_sizes_file_name(const std::filesystem::path &output_prefix) {
  assert(!output_prefix.empty());
  return fmt::format("{}.chrom.sizes", output_prefix);
}

std::filesystem::path generate_output_file_name(const std::filesystem::path &output_prefix,
                                                const hictk::Chromosome &chrom1,
                                                const hictk::Chromosome &chrom2) {
  return fmt::format("{}.{}.{}.parquet", output_prefix.string(), chrom1.name(), chrom2.name());
}

}  // namespace nchg
