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

#include <blockingconcurrentqueue.h>
#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <hictk/file.hpp>
#include <hictk/fmt/pixel.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <variant>

#ifndef _WIN32
#include <csignal>
#endif

// clang-format off
// As of HighFive 2.9.0, these headers must be included after HighFive/hictk,
// otherwise this source file fails to compile with MSVC
#include <boost/process/child.hpp>
#include <boost/process/pipe.hpp>
#include <boost/process/io.hpp>
// clang-format on

#include "nchg/common.hpp"
#include "nchg/config.hpp"
#include "nchg/nchg.hpp"
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
                 "observed_count\t"
                 "expected_count\t"
                 "odds_ratio\t"
                 "omega\n"));
  fflush(stdout);
}

[[nodiscard]] static std::string_view truncate_bed3_record(std::string_view record,
                                                           char sep = '\t') {
  const auto pos1 = record.find(sep);
  if (pos1 == std::string_view::npos) {
    throw std::runtime_error("invalid bed record, expected 3 tokens, found 1");
  }
  const auto pos2 = record.find('\t', pos1 + 1);
  if (pos2 == std::string_view::npos) {
    throw std::runtime_error("invalid bed record, expected 3 tokens, found 2");
  }
  const auto pos3 = record.find('\t', pos2 + 1);

  return record.substr(0, pos3);
}

[[nodiscard]] static std::vector<hictk::GenomicInterval> parse_domains(
    const hictk::Reference &chroms, const std::filesystem::path &path, std::string_view chrom1,
    std::string_view chrom2) {
  SPDLOG_INFO(FMT_STRING("reading domains from {}..."), path);
  std::vector<hictk::GenomicInterval> domains{};
  std::string buffer{};

  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  try {
    fs.open(path);

    for (std::size_t i = 1; std::getline(fs, buffer); ++i) {
      if (buffer.empty()) {
        continue;
      }

      if (buffer.back() == '\r') {
        buffer.resize(buffer.size() - 1);
      }

      try {
        const auto record = truncate_bed3_record(buffer);
        auto domain = hictk::GenomicInterval::parse_bed(chroms, record);

        if (chrom1 != "all") {
          assert(chrom2 != "all");
          if (domain.chrom().name() != chrom1 && domain.chrom().name() != chrom2) {
            continue;
          }
        }

        domains.emplace_back(std::move(domain));
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("found an invalid record at line {} of file {}: {}"), i, path, e.what()));
      }
    }

  } catch (const std::exception &e) {
    if (!fs.eof()) {
      throw;
    }
  }

  std::sort(domains.begin(), domains.end());
  SPDLOG_INFO(FMT_STRING("read {} domains from {}..."), domains.size(), path);
  return domains;
}

template <typename FilePtr, typename File = remove_cvref_t<decltype(*std::declval<FilePtr>())>>
[[nodiscard]] static NCHG<File> init_nchg(const FilePtr &f, const ComputePvalConfig &c) {
  assert(c.chrom1 != "all");
  assert(c.chrom2 != "all");
  assert(!c.cis_only);
  assert(!c.trans_only);

  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  if (!c.path_to_expected_values.empty()) {
    SPDLOG_INFO(FMT_STRING("reading expected values from {}..."), c.path_to_expected_values);
    return NCHG(f, chrom1, chrom2, ExpectedValues<File>::deserialize(c.path_to_expected_values));
  }

  return NCHG<File>(f, chrom1, chrom2, c.mad_max, c.min_delta, c.max_delta);
}

template <typename FilePtr>
[[nodiscard]] static std::size_t process_domains(const FilePtr &f, const ComputePvalConfig &c) {
  assert(std::filesystem::exists(c.path_to_domains));

  const auto domains = parse_domains(f->chromosomes(), c.path_to_domains, c.chrom1, c.chrom2);

  if (domains.empty()) {
    return 0;
  }

  const auto nchg = init_nchg(f, c);

  if (c.write_header) {
    print_header();
  }

  std::size_t num_records = 0;
  for (std::size_t i = 0; i < domains.size(); ++i) {
    for (std::size_t j = i; j < domains.size(); ++j) {
      const auto &d1 = domains[i];
      const auto &d2 = domains[j];

      if (c.chrom1 != "all" && (d1.chrom() != c.chrom1 || d2.chrom() != c.chrom2)) {
        continue;
      }

      const auto s = nchg.compute(d1, d2, c.bad_bin_fraction);
      fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval, s.pixel.count,
                 s.expected, s.odds_ratio, s.omega);
      ++num_records;
    }
  }
  return num_records;
}

template <typename FilePtr>
[[nodiscard]] static std::size_t process_one_chromosome_pair(const FilePtr &f,
                                                             const ComputePvalConfig &c) {
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);
  auto nchg = init_nchg(f, c);

  if (c.write_header) {
    print_header();
  }

  std::size_t num_records = 0;
  std::for_each(nchg.begin(chrom1, chrom2), nchg.end(chrom1, chrom2), [&](const auto &s) {
    ++num_records;
    fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval, s.pixel.count,
               s.expected, s.odds_ratio, s.omega);
  });
  return num_records;
}

[[nodiscard]] static std::size_t run_nchg_compute_worker(const ComputePvalConfig &c) {
  assert(c.chrom1 != "all");
  // clang-format off
  using FilePtr =
      std::variant<
          std::shared_ptr<const hictk::cooler::File>,
          std::shared_ptr<const hictk::hic::File>>;
  // clang-format on

  const auto f = [&]() -> FilePtr {
    hictk::File f_(c.path_to_hic.string(), c.resolution);
    return {std::visit(
        [&](auto &&ff) {
          using FileT = std::remove_reference_t<decltype(ff)>;
          return FilePtr{std::make_shared<const FileT>(std::forward<FileT>(ff))};
        },
        f_.get())};
  }();

  return std::visit(
      [&](const auto &f_) -> std::size_t {
        if (!c.path_to_domains.empty()) {
          return process_domains(f_, c);
        }
        return process_one_chromosome_pair(f_, c);
      },
      f);
}

[[nodiscard]] static std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>>
init_cis_chromosomes(const hictk::Reference &chroms) {
  std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> buffer{};

  for (const auto &chrom : chroms) {
    if (chrom.is_all()) {
      continue;
    }
    buffer.emplace_back(chrom, chrom);
  }

  return buffer;
}

[[nodiscard]] static std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>>
init_trans_chromosomes(const hictk::Reference &chroms) {
  std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> buffer{};

  for (const auto &chrom1 : chroms) {
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1.id() + 1; chrom2_id < chroms.size(); ++chrom2_id) {
      buffer.emplace_back(chrom1, chroms.at(chrom2_id));
    }
  }
  return buffer;
}

static void io_worker(moodycamel::BlockingConcurrentQueue<std::string> &msg_queue,
                      const std::atomic<bool> &early_return,
                      const std::atomic<std::size_t> &proc_completed,
                      const std::atomic<std::size_t> &proc_submitted,
                      const std::atomic<std::size_t> &msg_submitted,
                      std::atomic<std::size_t> &msg_received) {
  try {
    SPDLOG_DEBUG("spawning IO thread");
    moodycamel::ConsumerToken ctok(msg_queue);
    std::string buffer{};
    while (!early_return && (proc_completed != proc_submitted || msg_submitted != msg_received)) {
      SPDLOG_DEBUG(FMT_STRING("[IO] reading message..."));
      const auto msg_dequeued =
          msg_queue.wait_dequeue_timed(ctok, buffer, std::chrono::milliseconds(10));
      if (msg_dequeued) {
        SPDLOG_DEBUG(FMT_STRING("[IO] message read successfully!"));
        ++msg_received;
        fmt::print(FMT_COMPILE("{}\n"), buffer);
      } else {
        SPDLOG_DEBUG(FMT_STRING("[IO] unable to read message, trying again in 10 ms"));
      }
    }
    SPDLOG_DEBUG("[IO] returning");
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred in the IO thread: {}"), e.what()));
  }
}

[[nodiscard]] static auto spawn_compute_process(const ComputePvalConfig &c,
                                                const std::string &chrom1,
                                                const std::string &chrom2) {
  SPDLOG_INFO(FMT_STRING("begin processing {}:{}"), chrom1, chrom2);

  std::vector<std::string> args{"compute",
                                "--chrom1",
                                chrom1,
                                "--chrom2",
                                chrom2,
                                "--no-write-header",
                                "--write-eof",
                                "--threads",
                                "1",
                                "--verbosity",
                                "2",
                                "--min-delta",
                                fmt::to_string(c.min_delta),
                                "--max-delta",
                                fmt::to_string(c.max_delta),
                                "--resolution",
                                fmt::to_string(c.resolution),
                                c.path_to_hic.string()};

  if (!c.path_to_domains.empty()) {
    args.emplace_back("--domains");
    args.emplace_back(c.path_to_domains.string());
  }

  if (!c.path_to_expected_values.empty()) {
    args.emplace_back("--expected-values");
    args.emplace_back(c.path_to_expected_values.string());
  }

  boost::process::ipstream pipe;
  boost::process::child proc(
      c.exec.string(), args,
      boost::process::std_in<boost::process::null, boost::process::std_out> pipe);
  SPDLOG_DEBUG(FMT_STRING("spawned worker process {}..."), proc.id());
  if (!proc.running()) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to spawn worker process: {} {}"),
                                         c.exec.string(), fmt::join(args, " ")));
  }

  return std::make_pair(std::move(pipe), std::move(proc));
}

static void consume_compute_process_output(
    boost::process::child &proc, boost::process::ipstream &pipe,
    [[maybe_unused]] std::string_view chrom1, [[maybe_unused]] std::string_view chrom2,
    moodycamel::BlockingConcurrentQueue<std::string> &msg_queue, std::atomic<bool> &early_return,
    std::atomic<std::size_t> &msg_submitted) {
  std::string line;

  std::size_t records_processed{};
  while (true) {
    std::getline(pipe, line);

    if (line == "__EOF__" || line == "__EOF__\r") {
      break;
    }

    SPDLOG_DEBUG(FMT_STRING("[{}] sending message..."), proc.id());
    while (!early_return && !msg_queue.try_enqueue(line)) {
      SPDLOG_DEBUG(FMT_STRING("[{}] sending message failed! Retrying in 10 ms..."), proc.id());
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    ++msg_submitted;
    ++records_processed;

    if (early_return) {
      SPDLOG_DEBUG(FMT_STRING("[{}] terminating process..."), proc.id());
      proc.terminate();
      return;
    }
  }
  SPDLOG_DEBUG(FMT_STRING("[{}] done reading output from process!"), proc.id());
  SPDLOG_INFO(FMT_STRING("done processing {}:{} ({} records)!"), chrom1, chrom2, records_processed);
}

template <typename PidT>
[[nodiscard]] static std::size_t register_process(
    [[maybe_unused]] boost::process::child &proc, [[maybe_unused]] std::atomic<PidT *> &pids,
    [[maybe_unused]] const std::atomic<std::size_t> &num_pids,
    [[maybe_unused]] std::mutex &pids_mtx) {
#ifdef _WIN32
  return 0;
#else
  const std::scoped_lock lck(pids_mtx);
  for (std::size_t i = 0; i < num_pids; ++i) {
    if (pids[i] == -1) {
      pids[i] = conditional_static_cast<PidT>(proc.id());
      return i;
    }
  }
  proc.terminate();
  throw std::runtime_error(fmt::format(FMT_STRING("unable to register process {}"), proc.id()));
#endif
}

template <typename PidT>
static void deregister_process([[maybe_unused]] std::size_t slot,
                               [[maybe_unused]] std::atomic<PidT *> &pids) {
#ifndef _WIN32
  pids[slot] = static_cast<PidT>(-1);
#endif
}

template <typename PidT>
static std::size_t process_queries(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const ComputePvalConfig &c, std::atomic<PidT *> &pids,
    const std::atomic<std::size_t> &num_pids) {
  BS::thread_pool tpool(conditional_static_cast<BS::concurrency_t>(c.threads + 1));
  moodycamel::BlockingConcurrentQueue<std::string> msg_queue(c.threads * 1'000);

  std::mutex pids_mtx{};

  std::atomic<std::size_t> msg_submitted{};
  std::atomic<std::size_t> msg_received{};

  const std::atomic<std::size_t> proc_submitted{chrom_pairs.size()};
  std::atomic<std::size_t> proc_completed{};
  std::atomic<bool> early_return{false};

  auto io = tpool.submit_task([&]() {
    io_worker(msg_queue, early_return, proc_completed, proc_submitted, msg_submitted, msg_received);
  });

  auto workers = tpool.submit_loop<std::size_t>(
      0, chrom_pairs.size(),
      [&](std::size_t i) {
        if (early_return) {
          return;
        }

        const auto chrom1 = std::string{chrom_pairs[i].first.name()};
        const auto chrom2 = std::string{chrom_pairs[i].second.name()};

        try {
          auto [pipe, proc] = spawn_compute_process(c, chrom1, chrom2);
          const auto pid_offset = register_process(proc, pids, num_pids, pids_mtx);

          consume_compute_process_output(proc, pipe, chrom1, chrom2, msg_queue, early_return,
                                         msg_submitted);

          proc.wait();
          deregister_process(pid_offset, pids);
          ++proc_completed;

          if (proc.exit_code() != 0) {
            early_return = true;
            throw std::runtime_error(
                fmt::format(FMT_STRING("child process terminated with code {}"), proc.exit_code()));
          }
        } catch (const std::exception &e) {
          early_return = true;
          throw std::runtime_error(
              fmt::format(FMT_STRING("error in the worker thread processing {}:{}: {}"), chrom1,
                          chrom2, e.what()));
        } catch (...) {
          early_return = true;
        }
      },
      chrom_pairs.size());

  workers.get();
  io.get();

  return msg_received.load();
}

template <typename PidT>
int run_nchg_compute(const ComputePvalConfig &c, std::atomic<PidT *> &pids,
                     const std::atomic<std::size_t> &num_pids) {
  const auto t0 = std::chrono::system_clock::now();

  if (c.chrom1 != "all") {
    assert(c.chrom2 != "all");
    const auto interactions_processed = run_nchg_compute_worker(c);
    if (c.write_eof_signal) {
      fmt::print(FMT_STRING("__EOF__\n"));
    }
    const auto t1 = std::chrono::system_clock::now();
    const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    SPDLOG_INFO(FMT_STRING("Processed {} interactions in {}s!"), interactions_processed,
                static_cast<double>(delta) / 1000.0);
    return 0;
  }

  const hictk::File f(c.path_to_hic.string(), c.resolution);
  std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> chrom_pairs{};
  if (c.cis_only) {
    chrom_pairs = init_cis_chromosomes(f.chromosomes());
  }
  if (c.trans_only) {
    chrom_pairs = init_trans_chromosomes(f.chromosomes());
  }

  if (c.chrom1 == "all" && !c.cis_only && !c.trans_only) {
    chrom_pairs = init_cis_chromosomes(f.chromosomes());
    const auto chrom_pairs2 = init_trans_chromosomes(f.chromosomes());
    std::copy(chrom_pairs2.begin(), chrom_pairs2.end(), std::back_inserter(chrom_pairs));
  }

  if (c.write_header) {
    print_header();
  }
  const auto interactions_processed = process_queries(chrom_pairs, c, pids, num_pids);

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  SPDLOG_INFO(FMT_STRING("Processed {} interactions in {}s!"), interactions_processed,
              static_cast<double>(delta) / 1000.0);

  return 0;
}

#ifdef _WIN32
int run_nchg_compute(const ComputePvalConfig &c, std::atomic<std::uint32_t *> &pids,
                     const std::atomic<std::size_t> &num_pids) {
  return run_nchg_compute<std::uint32_t>(c, pids, num_pids);
}
#else
int run_nchg_compute(const ComputePvalConfig &c, std::atomic<pid_t *> &pids,
                     const std::atomic<std::size_t> &num_pids) {
  return run_nchg_compute<pid_t>(c, pids, num_pids);
}
#endif
}  // namespace nchg
