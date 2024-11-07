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



#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/io/file.h>
#include <hictk/file.hpp>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/std.h>
#include <parquet/stream_reader.h>
#include <parquet/stream_writer.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <hictk/fmt/pixel.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <type_traits>
#include <variant>

// clang-format off
// As of HighFive 2.9.0, these headers must be included after HighFive/hictk,
// otherwise this source file fails to compile with MSVC
#include <boost/process/child.hpp>
#include <boost/process/pipe.hpp>
#include <boost/process/io.hpp>
// clang-format on

#include "nchg/common.hpp"
#include "nchg/concepts.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/io.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

[[nodiscard]] static std::vector<hictk::GenomicInterval> parse_domains(
    const hictk::Reference &chroms, const std::filesystem::path &path, std::string_view chrom1,
    std::string_view chrom2) {
  SPDLOG_INFO("[{}:{}] reading domains from {}...", chrom1, chrom2, path);
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
        throw std::runtime_error(
            fmt::format("found an invalid record at line {} of file {}: {}", i, path, e.what()));
      }
    }

  } catch (const std::exception &) {
    if (!fs.eof()) {
      throw;
    }
  }

  std::ranges::sort(domains);
  SPDLOG_INFO("[{}:{}] read {} domains from {}...", chrom1, chrom2, domains.size(), path);
  return domains;
}

[[nodiscard]] static NCHG init_nchg(const std::shared_ptr<const hictk::File> &f,
                                    const std::optional<ExpectedValues> &expected_values,
                                    const ComputePvalConfig &c) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(c.compute_cis);
  assert(c.compute_trans);

  const auto &chrom1 = f->chromosomes().at(*c.chrom1);
  const auto &chrom2 = f->chromosomes().at(*c.chrom2);

  if (expected_values.has_value()) {
    return {f, chrom1, chrom2, *expected_values};
  }

  if (!c.path_to_expected_values.empty()) {
    SPDLOG_INFO("reading expected values from {}...", c.path_to_expected_values);
    return {f, chrom1, chrom2, ExpectedValues::deserialize(c.path_to_expected_values)};
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

  return {f, chrom1, chrom2, evs};
}

static void write_chrom_sizes_to_file(const hictk::Reference &chroms,
                                      const std::filesystem::path &path, bool force) {
  try {
    if (force) {
      std::filesystem::remove(path);  // NOLINT
    }

    if (std::filesystem::exists(path)) {
      throw std::runtime_error(
          fmt::format("Refusing to overwrite file {}. Pass --force to overwrite.", path));
    }

    const auto output_dir = path.parent_path();
    if (!output_dir.empty() && !std::filesystem::exists(output_dir)) {
      std::filesystem::create_directories(output_dir);
    }

    std::ofstream fs{};
    fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);
    fs.open(path);

    for (const auto &chrom : chroms) {
      fmt::print(fs, "{}\t{}\n", chrom.name(), chrom.size());
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to write chromosomes to file {}: {}", path, e.what()));
  }
}

[[nodiscard]] static std::size_t process_domains(
    const std::shared_ptr<const hictk::File> &f,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  assert(std::filesystem::exists(c.path_to_domains));
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(!c.output_path.empty());

  SPDLOG_INFO("[{}:{}] begin processing domains from {}...", *c.chrom1, *c.chrom2,
              c.path_to_domains);

  const auto writer = init_parquet_file_writer<NCHGResult>(
      f->chromosomes(), c.output_path, c.force, c.compression_method, c.compression_lvl, c.threads);

  const auto domains = parse_domains(f->chromosomes(), c.path_to_domains, *c.chrom1, *c.chrom2);

  if (domains.empty()) {
    return 0;
  }

  const auto nchg = init_nchg(f, expected_values, c);

  constexpr std::size_t batch_size = 1'000'000;
  RecordBatchBuilder builder(f->bins().chromosomes());

  std::size_t num_records = 0;
  for (std::size_t i = 0; i < domains.size(); ++i) {
    for (std::size_t j = i; j < domains.size(); ++j) {
      const auto &d1 = domains[i];
      const auto &d2 = domains[j];

      if (c.chrom1.has_value() && (d1.chrom() != *c.chrom1 || d2.chrom() != *c.chrom2)) {
        continue;
      }

      const auto s = nchg.compute(d1, d2, c.bad_bin_fraction);

      if (builder.size() == batch_size) [[unlikely]] {
        builder.write(*writer);
      }

      if (std::isfinite(s.odds_ratio) && s.omega != 0) [[likely]] {
        builder.append(s);
      }

      ++num_records;
    }
  }

  if (builder.size() != 0) {
    builder.write(*writer);
  }
  return num_records;
}

[[nodiscard]] static std::size_t process_one_chromosome_pair(
    const std::shared_ptr<const hictk::File> &f,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(!c.output_path.empty());

  SPDLOG_INFO("[{}:{}] begin processing interactions...", *c.chrom1, *c.chrom2);

  const auto &chrom1 = f->chromosomes().at(*c.chrom1);
  const auto &chrom2 = f->chromosomes().at(*c.chrom2);
  const auto nchg = init_nchg(f, expected_values, c);

  const auto writer = init_parquet_file_writer<NCHGResult>(
      f->chromosomes(), c.output_path, c.force, c.compression_method, c.compression_lvl, c.threads);

  constexpr std::size_t batch_size = 1'000'000;
  RecordBatchBuilder builder{f->chromosomes()};

  std::size_t num_records = 0;
  auto first = nchg.begin(chrom1, chrom2);
  auto last = nchg.end(chrom1, chrom2);

  std::visit(
      [&]<typename It>(It &it1) {
        auto &it2 = std::get<It>(last);
        std::for_each(it1, it2, [&](const auto &s) {
          ++num_records;

          if (builder.size() == batch_size) [[unlikely]] {
            builder.write(*writer);
          }

          if (std::isfinite(s.odds_ratio) && s.omega != 0) [[likely]] {
            builder.append(s);
          }
        });
      },
      first);

  if (builder.size() != 0) {
    builder.write(*writer);
  }
  return num_records;
}

[[nodiscard]] static std::size_t run_nchg_compute_worker(
    const ComputePvalConfig &c, const std::optional<ExpectedValues> &expected_values = {}) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());

  const auto f = std::make_shared<const hictk::File>(c.path_to_hic.string(), c.resolution);

  if (!c.path_to_domains.empty()) {
    return process_domains(f, expected_values, c);
  }

  return process_one_chromosome_pair(f, expected_values, c);
}

using ChromosomePairs = std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>>;

[[nodiscard]] static ChromosomePairs init_cis_chromosomes(const hictk::Reference &chroms) {
  ChromosomePairs buffer{};

  for (const auto &chrom : chroms) {
    if (chrom.is_all()) [[unlikely]] {
      continue;
    }
    buffer.emplace_back(chrom, chrom);
  }

  return buffer;
}

[[nodiscard]] static ChromosomePairs init_trans_chromosomes(const hictk::Reference &chroms) {
  ChromosomePairs buffer{};

  for (const auto &chrom1 : chroms) {
    if (chrom1.is_all()) [[unlikely]] {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1.id() + 1; chrom2_id < chroms.size(); ++chrom2_id) {
      buffer.emplace_back(chrom1, chroms.at(chrom2_id));
    }
  }
  return buffer;
}

[[nodiscard]] static boost::process::child spawn_compute_process(const ComputePvalConfig &c,
                                                                 const hictk::Chromosome &chrom1,
                                                                 const hictk::Chromosome &chrom2) {
  assert(c.output_prefix.empty());
  assert(!c.output_path.empty());

  SPDLOG_INFO("[{}:{}]: begin processing...", chrom1.name(), chrom2.name());

  std::vector<std::string> args{
      "compute",
      "--chrom1",
      std::string{chrom1.name()},
      "--chrom2",
      std::string{chrom2.name()},
      "--threads",
      "1",
      "--resolution",
      fmt::to_string(c.resolution),
      "--bad-bin-fraction",
      fmt::to_string(c.bad_bin_fraction),
      "--compression-level",
      fmt::to_string(c.compression_lvl),
      "--compression-method",
      c.compression_method,
      "--verbosity",
      "2",
      c.path_to_hic.string(),
      c.output_path.string(),
  };

  if (!c.path_to_domains.empty()) {
    args.emplace_back("--domains");
    args.emplace_back(c.path_to_domains.string());
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

  for (std::size_t attempt = 0; attempt < 10; ++attempt) {
    boost::process::child proc(
        c.exec.string(), args,
        boost::process::std_in<boost::process::null, boost::process::std_out> boost::process::null);
    if (proc.running() || proc.exit_code() == 0) {
      SPDLOG_DEBUG("[{}:{}]: spawned worker process {}...", chrom1.name(), chrom2.name(),
                   proc.id());
      return proc;
    }
    SPDLOG_WARN("[{}:{}]: spawning worker process {} failed (attempt {}/10)...", chrom1.name(),
                chrom2.name(), proc.id(), attempt + 1);
    proc.terminate();
  }

  throw std::runtime_error(
      fmt::format("failed to spawn worker process: {} {}", c.exec.string(), fmt::join(args, " ")));
}

[[nodiscard]] static std::filesystem::path generate_chrom_sizes_file_name(
    const std::filesystem::path &output_prefix) {
  assert(!output_prefix.empty());
  return fmt::format("{}.chrom.sizes", output_prefix);
}

[[nodiscard]] static std::filesystem::path generate_output_file_name(
    const std::filesystem::path &output_prefix, const hictk::Chromosome &chrom1,
    const hictk::Chromosome &chrom2) {
  return fmt::format("{}.{}.{}.parquet", output_prefix.string(), chrom1.name(), chrom2.name());
}

[[nodiscard]] static std::size_t worker_fx(const hictk::Chromosome &chrom1,
                                           const hictk::Chromosome &chrom2,
                                           const ComputePvalConfig &config,
                                           bool trans_expected_values_avail,
                                           std::atomic<bool> &early_return) {
  if (early_return) {
    return 0;
  }

  try {
    auto child_config = config;

    child_config.output_path =
        generate_output_file_name(config.output_prefix.string(), chrom1, chrom2);
    child_config.output_prefix.clear();

    if (!trans_expected_values_avail && chrom1 != chrom2) {
      child_config.path_to_expected_values.clear();
    }

    auto proc = spawn_compute_process(child_config, chrom1, chrom2);
    proc.wait();

    if (proc.exit_code() != 0) {
      early_return = true;
      throw std::runtime_error(
          fmt::format("child process terminated with code {}", proc.exit_code()));
    }

    std::shared_ptr<arrow::io::ReadableFile> fp;
    PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(child_config.output_path));
    const auto records_processed =
        parquet::StreamReader{parquet::ParquetFileReader::Open(fp)}.num_rows();

    SPDLOG_INFO("done processing {}:{} ({} records)!", chrom1.name(), chrom2.name(),
                records_processed);
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
    const ComputePvalConfig &c, const std::optional<ExpectedValues> &expected_values) {
  auto base_config = c;
  if (c.path_to_expected_values.empty() && expected_values.has_value()) {
    assert(!base_config.output_prefix.empty());
    const auto tmpdir = base_config.output_prefix.parent_path() / "tmp/";
    if (std::filesystem::exists(tmpdir) && !c.force) {
      throw std::runtime_error(
          fmt::format("refusing to overwrite existing temporary folder: \"{}\". "
                      "Pass --force to overwrite.",
                      tmpdir));
    }

    base_config.path_to_expected_values =
        tmpdir / fmt::format("{}_expected_values.h5", base_config.output_prefix.stem().string());

    if (!base_config.force && std::filesystem::exists(base_config.path_to_expected_values)) {
      throw std::runtime_error(
          fmt::format("Refusing to overwrite file {}. Pass --force to overwrite.",
                      base_config.path_to_expected_values));
    }

    std::filesystem::create_directories(tmpdir);
    std::filesystem::remove(base_config.path_to_expected_values);  // NOLINT
    expected_values->serialize(base_config.path_to_expected_values);
  }

  return base_config;
}

static std::size_t process_queries_mt(
    BS::thread_pool &tpool,
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  std::atomic<bool> early_return{false};

  const auto user_provided_expected_values = !c.path_to_expected_values.empty();
  const auto base_config = init_base_config(c, expected_values);

  BS::multi_future<std::size_t> workers(chrom_pairs.size());
  for (std::size_t i = 0; i < workers.size(); ++i) {
    workers[i] = tpool.submit_task([&, i] {
      const auto &[chrom1, chrom2] = chrom_pairs[i];
      SPDLOG_DEBUG("submitting task {}/{} ({}:{})...", i + 1, workers.size(), chrom1.name(),
                   chrom2.name());
      return worker_fx(chrom1, chrom2, base_config, user_provided_expected_values, early_return);
    });
  }

  std::size_t num_records = 0;
  for (const auto result : workers.get()) {
    num_records += result;
  }

  return num_records;
}

static std::size_t process_queries_st(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
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
                                   ? run_nchg_compute_worker(config, expected_values)
                                   : run_nchg_compute_worker(config);
      tot_num_records += num_records;

      const auto t1 = std::chrono::steady_clock::now();
      SPDLOG_INFO("[{}:{}] processed {} records in {}!", chrom1.name(), chrom2.name(), num_records,
                  format_duration(t1 - t0));
      SPDLOG_INFO("[{}:{}] {} records have been written to file \"{}\"", chrom1.name(),
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

static std::optional<ExpectedValues> init_cis_expected_values(const ComputePvalConfig &c) {
  assert(c.compute_cis || c.chrom1 == c.chrom2);

  SPDLOG_INFO("initializing expected values for cis matrices...");
  const auto f = std::make_shared<hictk::File>(c.path_to_hic.string(), c.resolution);

  const auto bin_mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  return {ExpectedValues::cis_only(
      f,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      bin_mask)};
}

static std::size_t process_queries(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  write_chrom_sizes_to_file(hictk::File(c.path_to_hic, c.resolution).chromosomes(),
                            generate_chrom_sizes_file_name(c.output_prefix), c.force);

  if (c.threads > 1) {
    auto num_threads = conditional_static_cast<BS::concurrency_t>(c.threads);
    if (c.threads > chrom_pairs.size()) {
      num_threads = conditional_static_cast<BS::concurrency_t>(chrom_pairs.size());
      SPDLOG_WARN(
          "number of threads specified through --threads exceeds the number of chromosome pairs to "
          "be processed: limiting concurrency to {} thread(s)",
          num_threads);
    }
    BS::thread_pool tpool(num_threads);
    return process_queries_mt(tpool, chrom_pairs, expected_values, c);
  }
  return process_queries_st(chrom_pairs, expected_values, c);
}

static void process_file_collisions(
    const std::filesystem::path &output_prefix,
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs, bool force) {
  std::vector<std::string> collisions{};
  std::size_t num_collisions = 0;

  constexpr std::size_t max_collisions_reported = 10;

  if (const auto chrom_sizes = generate_chrom_sizes_file_name(output_prefix); force) {
    const auto removed = std::filesystem::remove(chrom_sizes);  // NOLINT
    if (removed) {
      SPDLOG_DEBUG("file \"{}\" has been deleted", chrom_sizes.string());
    }
  } else if (std::filesystem::exists(chrom_sizes)) {
    collisions.emplace_back(chrom_sizes.string());
    ++num_collisions;
  }

  for (const auto &[chrom1, chrom2] : chrom_pairs) {
    const auto output_file_name = generate_output_file_name(output_prefix, chrom1, chrom2);

    if (force) {
      // We are removing files eagerly to avoid scenarios where:
      // - the user is running NCHG compute on a dirty folder using --force
      // - NCHG crashes or exits with an error message
      // - the user does not notice the error, but observes that all output files are there and
      //   assumes that all is good

      const auto removed = std::filesystem::remove(output_file_name);  // NOLINT
      if (removed) {
        SPDLOG_DEBUG("file \"{}\" has been deleted", output_file_name.string());
      }
    } else if (std::filesystem::exists(output_file_name)) {
      ++num_collisions;
      if (collisions.size() < max_collisions_reported) {
        collisions.emplace_back(output_file_name.string());
      }
    }
  }

  if (num_collisions != 0) {
    assert(!collisions.empty());
    const std::string suffix =
        num_collisions > collisions.size()
            ? fmt::format("\n - and {} more", num_collisions - collisions.size())
            : "";

    throw std::runtime_error(fmt::format(
        "refusing to overwrite the following file(s). Pass --force to overwrite.\n - {}{}",
        fmt::join(collisions, "\n - "), suffix));
  }
}

static void validate_expected_values(
    const ExpectedValues &expected_values, const std::filesystem::path &path_to_expected_values,
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    std::uint32_t resolution) {
  if (expected_values.resolution() != resolution) {
    throw std::runtime_error(
        fmt::format("mismatch in file resolution: expected values have been computed "
                    "for {}bp resolution but given Hi-C matrix has {}bp resolution",
                    expected_values.resolution(), resolution));
  }

  std::vector<std::string> missing_values{};
  std::size_t num_missing_values = 0;

  constexpr std::size_t max_missing_values_reported = 10;

  for (const auto &[chrom1, chrom2] : chrom_pairs) {
    try {
      if (chrom1 == chrom2) {
        std::ignore = expected_values.expected_values(chrom1);
      } else if (!path_to_expected_values.empty()) {
        std::ignore = expected_values.expected_value(chrom1, chrom2);
      }
    } catch (const std::out_of_range &) {
      ++num_missing_values;
      if (missing_values.size() < max_missing_values_reported) {
        missing_values.emplace_back(fmt::format("{}:{}", chrom1.name(), chrom2.name()));
      }
    }
  }

  if (num_missing_values != 0) {
    assert(!missing_values.empty());
    const std::string suffix =
        num_missing_values > missing_values.size()
            ? fmt::format("\n - and {} more", num_missing_values - missing_values.size())
            : "";

    throw std::runtime_error(
        fmt::format("user-specified expected value file (\"{}\") does not contain information for "
                    "the following chromosome pair(s):\n - {}{}",
                    path_to_expected_values, fmt::join(missing_values, "\n - "), suffix));
  }
}

[[nodiscard]] static auto generate_execution_plan(const ComputePvalConfig &c) {
  struct Plan {
    std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> chrom_pairs{};
    std::optional<ExpectedValues> expected_values{};
  };

  Plan plan{};

  const hictk::File f(c.path_to_hic.string(), c.resolution);

  if (f.bins().type() != hictk::BinTable::Type::fixed) {
    throw std::runtime_error("only file with uniform bin sizes are supported.");
  }

  if (c.compute_cis) {
    plan.chrom_pairs = init_cis_chromosomes(f.chromosomes());
  }

  if (c.compute_trans) {
    std::ranges::copy(init_trans_chromosomes(f.chromosomes()),
                      std::back_inserter(plan.chrom_pairs));
  }

  process_file_collisions(c.output_prefix, plan.chrom_pairs, c.force);

  if (!c.path_to_expected_values.empty()) {
    assert(!plan.expected_values.has_value());
    plan.expected_values = ExpectedValues::deserialize(c.path_to_expected_values);
  } else if (c.compute_cis) {
    plan.expected_values = init_cis_expected_values(c);
  }

  if (plan.expected_values.has_value()) {
    validate_expected_values(*plan.expected_values, c.path_to_expected_values, plan.chrom_pairs,
                             f.resolution());
  }

  return plan;
}

int run_nchg_compute(const ComputePvalConfig &c) {
  try {
    const auto t0 = std::chrono::steady_clock::now();

    if (c.chrom1.has_value()) {
      assert(c.chrom2.has_value());
      assert(!c.output_path.empty());
      const auto interactions_processed = run_nchg_compute_worker(c);
      const auto t1 = std::chrono::steady_clock::now();
      SPDLOG_INFO("[{}:{}] processed {} records in {}!", *c.chrom1, *c.chrom2,
                  interactions_processed, format_duration(t1 - t0));

      SPDLOG_INFO("records for {}:{} have been written to file \"{}\"", *c.chrom1, *c.chrom2,
                  c.output_path.string());
      return 0;
    }

    const auto &[chrom_pairs, expected_values] = generate_execution_plan(c);

    const auto interactions_processed = process_queries(chrom_pairs, expected_values, c);
    std::filesystem::remove_all(c.tmpdir);  // NOLINT

    const auto t1 = std::chrono::steady_clock::now();
    if (interactions_processed == 0) {
      SPDLOG_WARN("no records have been processed. Is this intended?");
    } else {
      SPDLOG_INFO("processed {} records in {}!", interactions_processed, format_duration(t1 - t0));
    }

    SPDLOG_INFO("{} file(s) have been created under prefix \"{}\"", chrom_pairs.size() + 1,
                c.output_prefix);

    return 0;
  } catch (...) {
    std::filesystem::remove_all(c.tmpdir);  // NOLINT
    throw;
  }
}

}  // namespace nchg
