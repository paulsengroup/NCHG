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
#include "nchg/config.hpp"
#include "nchg/io.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {

[[nodiscard]] static std::vector<hictk::GenomicInterval> parse_domains(
    const hictk::Reference &chroms, const std::filesystem::path &path, std::string_view chrom1,
    std::string_view chrom2) {
  SPDLOG_INFO(FMT_STRING("[{}:{}] reading domains from {}..."), chrom1, chrom2, path);
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

  } catch (const std::exception &) {
    if (!fs.eof()) {
      throw;
    }
  }

  std::sort(domains.begin(), domains.end());
  SPDLOG_INFO(FMT_STRING("[{}:{}] read {} domains from {}..."), chrom1, chrom2, domains.size(),
              path);
  return domains;
}

template <typename File>
  requires HictkSingleResFile<File>
[[nodiscard]] static NCHG<File> init_nchg(
    const std::shared_ptr<const File> &f,
    const std::optional<ExpectedValues<File>> &expected_values, const ComputePvalConfig &c) {
  assert(c.chrom1 != "all");
  assert(c.chrom2 != "all");
  assert(!c.cis_only);
  assert(!c.trans_only);

  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  if (expected_values.has_value()) {
    return NCHG(f, chrom1, chrom2, *expected_values);
  }

  if (!c.path_to_expected_values.empty()) {
    SPDLOG_INFO(FMT_STRING("reading expected values from {}..."), c.path_to_expected_values);
    return NCHG(f, chrom1, chrom2, ExpectedValues<File>::deserialize(c.path_to_expected_values));
  }

  const auto bin_mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  const auto evs = ExpectedValues<File>::chromosome_pair(
      f, chrom1, chrom2,
      {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
       c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
       c.interpolation_qtile, c.interpolation_window_size},
      bin_mask);

  return NCHG(f, chrom1, chrom2, evs);
}

static void write_chrom_sizes_to_file(const hictk::Reference &chroms,
                                      const std::filesystem::path &path, bool force) {
  try {
    if (force) {
      std::filesystem::remove(path);
    }

    if (std::filesystem::exists(path)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), path));
    }

    std::ofstream fs{};
    fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);
    fs.open(path);

    for (const auto &chrom : chroms) {
      fmt::print(fs, FMT_STRING("{}\t{}\n"), chrom.name(), chrom.size());
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to write chromosomes to file {}: {}"), path, e.what()));
  }
}

template <typename File>
  requires HictkSingleResFile<File>
[[nodiscard]] static std::size_t process_domains(
    const std::shared_ptr<const File> &f,
    const std::optional<ExpectedValues<File>> &expected_values, const ComputePvalConfig &c) {
  assert(std::filesystem::exists(c.path_to_domains));

  SPDLOG_INFO(FMT_STRING("[{}:{}] begin processing domains from {}..."), c.chrom1, c.chrom2,
              c.path_to_domains);

  auto writer =
      init_parquet_file_writer<NCHGResult>(f->chromosomes(), c.output_prefix, c.force,
                                           c.compression_method, c.compression_lvl, c.threads);

  const auto domains = parse_domains(f->chromosomes(), c.path_to_domains, c.chrom1, c.chrom2);

  if (domains.empty()) {
    return 0;
  }

  const auto nchg = init_nchg(f, expected_values, c);

  std::size_t batch_size = 1'000'000;
  RecordBatchBuilder builder(f->bins().chromosomes());

  std::size_t num_records = 0;
  for (std::size_t i = 0; i < domains.size(); ++i) {
    for (std::size_t j = i; j < domains.size(); ++j) {
      const auto &d1 = domains[i];
      const auto &d2 = domains[j];

      if (c.chrom1 != "all" && (d1.chrom() != c.chrom1 || d2.chrom() != c.chrom2)) {
        continue;
      }

      const auto s = nchg.compute(d1, d2, c.bad_bin_fraction);

      if (builder.size() == batch_size) {
        builder.write(*writer);
      }

      if (std::isfinite(s.odds_ratio) && s.omega != 0) {
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

template <typename File>
  requires HictkSingleResFile<File>
[[nodiscard]] static std::size_t process_one_chromosome_pair(
    const std::shared_ptr<const File> &f,
    const std::optional<ExpectedValues<File>> &expected_values, const ComputePvalConfig &c) {
  SPDLOG_INFO(FMT_STRING("[{}:{}] begin processing interactions..."), c.chrom1, c.chrom2);

  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);
  auto nchg = init_nchg(f, expected_values, c);

  auto writer =
      init_parquet_file_writer<NCHGResult>(f->chromosomes(), c.output_prefix, c.force,
                                           c.compression_method, c.compression_lvl, c.threads);

  const std::size_t batch_size = 1'000'000;
  RecordBatchBuilder builder{f->chromosomes()};

  std::size_t num_records = 0;
  std::for_each(nchg.begin(chrom1, chrom2), nchg.end(chrom1, chrom2), [&](const auto &s) {
    ++num_records;

    if (builder.size() == batch_size) {
      builder.write(*writer);
    }

    if (std::isfinite(s.odds_ratio) && s.omega != 0) {
      builder.append(s);
    }
  });

  if (builder.size() != 0) {
    builder.write(*writer);
  }
  return num_records;
}

// clang-format off
using HiCFilePtr =
    std::variant<
        std::shared_ptr<const hictk::cooler::File>,
        std::shared_ptr<const hictk::hic::File>>;
// clang-format on

[[nodiscard]] static HiCFilePtr open_file_ptr(const std::filesystem::path &path,
                                              std::uint32_t resolution) {
  hictk::File f_(path.string(), resolution);
  return {std::visit(
      [&](auto &&ff) {
        using FileT = std::remove_reference_t<decltype(ff)>;
        return HiCFilePtr{std::make_shared<const FileT>(std::forward<FileT>(ff))};
      },
      f_.get())};
}

template <typename File = hictk::cooler::File>
  requires HictkSingleResFile<File>
[[nodiscard]] static std::size_t run_nchg_compute_worker(
    const ComputePvalConfig &c, const std::optional<ExpectedValues<File>> &expected_values = {}) {
  assert(c.chrom1 != "all");

  const auto f = open_file_ptr(c.path_to_hic, c.resolution);

  return std::visit(
      [&]<typename FilePtr>(const FilePtr &f_) -> std::size_t {
        using UnderlyingFile = std::remove_cvref_t<typename FilePtr::element_type>;
        if constexpr (!std::is_same_v<File, UnderlyingFile>) {
          std::optional<ExpectedValues<UnderlyingFile>> expected_values_{};
          if (expected_values) {
            expected_values_ = expected_values->template cast<UnderlyingFile>();
          }

          return run_nchg_compute_worker(c, expected_values_);
        } else {
          if (!c.path_to_domains.empty()) {
            return process_domains(f_, expected_values, c);
          }

          return process_one_chromosome_pair(f_, expected_values, c);
        }
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

[[nodiscard]] static boost::process::child spawn_compute_process(
    const ComputePvalConfig &c, const std::filesystem::path &output_path,
    const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) {
  SPDLOG_INFO(FMT_STRING("begin processing {}:{}"), chrom1.name(), chrom2.name());

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
      fmt::to_string(c.compression_method),
      "--verbosity",
      "2",
      c.path_to_hic.string(),
      output_path.string(),
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
    if (proc.running()) {
      SPDLOG_DEBUG(FMT_STRING("spawned worker process {}..."), proc.id());
      return proc;
    }
    SPDLOG_DEBUG(FMT_STRING("spawning worker process failed (attempt {}/10)..."), proc.id(),
                 attempt + 1);
    proc.terminate();
  }

  throw std::runtime_error(fmt::format(FMT_STRING("failed to spawn worker process: {} {}"),
                                       c.exec.string(), fmt::join(args, " ")));
}

template <typename File>
  requires HictkSingleResFile<File>
static std::size_t process_queries_mt(
    BS::thread_pool &tpool,
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const std::optional<ExpectedValues<File>> &expected_values, const ComputePvalConfig &c) {
  std::atomic<bool> early_return{false};

  auto config = c;
  if (c.path_to_expected_values.empty() && expected_values.has_value()) {
    config.path_to_expected_values =
        fmt::format(FMT_STRING("{}_expected_values.h5"), config.output_prefix.string());

    if (!config.force && std::filesystem::exists(config.path_to_expected_values)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."),
                      config.path_to_expected_values));
    }

    std::filesystem::remove(config.path_to_expected_values);

    expected_values->serialize(config.path_to_expected_values);
  }

  auto workers = tpool.submit_blocks(
      std::size_t{0}, chrom_pairs.size(),
      [&](std::size_t i1, [[maybe_unused]] std::size_t i2) -> std::size_t {
        if (early_return) {
          return 0;
        }

        const auto &chrom1 = chrom_pairs[i1].first;
        const auto &chrom2 = chrom_pairs[i1].second;

        try {
          const auto output_path =
              fmt::format(FMT_STRING("{}.{}.{}.parquet"), config.output_prefix.string(),
                          chrom1.name(), chrom2.name());

          const auto trans_expected_values_avail = !c.path_to_expected_values.empty();

          auto child_config = config;
          if (!trans_expected_values_avail && chrom1 != chrom2) {
            child_config.path_to_expected_values.clear();
          }

          auto proc = spawn_compute_process(child_config, output_path, chrom1, chrom2);
          proc.wait();

          if (proc.exit_code() != 0) {
            early_return = true;
            throw std::runtime_error(
                fmt::format(FMT_STRING("child process terminated with code {}"), proc.exit_code()));
          }

          std::shared_ptr<arrow::io::ReadableFile> fp;
          PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(output_path));
          const auto records_processed =
              parquet::StreamReader{parquet::ParquetFileReader::Open(fp)}.num_rows();

          SPDLOG_INFO(FMT_STRING("done processing {}:{} ({} records)!"), chrom1.name(),
                      chrom2.name(), records_processed);
          return static_cast<std::size_t>(records_processed);

        } catch (const std::exception &e) {
          early_return = true;
          throw std::runtime_error(
              fmt::format(FMT_STRING("error in the worker thread processing {}:{}: {}"),
                          chrom1.name(), chrom2.name(), e.what()));
        } catch (...) {
          early_return = true;
          throw;
        }
      },
      chrom_pairs.size());

  std::size_t num_records = 0;
  const auto results = workers.get();
  for (const auto &res : results) {
    num_records += res;
  }

  return num_records;
}

template <typename File>
  requires HictkSingleResFile<File>
static std::size_t process_queries_st(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const std::optional<ExpectedValues<File>> &expected_values, const ComputePvalConfig &c) {
  std::size_t tot_num_records = 0;
  for (const auto &pair : chrom_pairs) {
    const auto &chrom1 = pair.first;
    const auto &chrom2 = pair.second;

    try {
      auto config = c;
      config.chrom1 = chrom1.name();
      config.chrom2 = chrom2.name();
      config.cis_only = false;
      config.trans_only = false;

      config.output_prefix = fmt::format(FMT_STRING("{}.{}.{}.parquet"), c.output_prefix.string(),
                                         chrom1.name(), chrom2.name());

      const auto t0 = std::chrono::system_clock::now();
      if (config.chrom1 == config.chrom2) {
        assert(expected_values.has_value());
      }
      const auto num_records = config.chrom1 == config.chrom2
                                   ? run_nchg_compute_worker(config, expected_values)
                                   : run_nchg_compute_worker(config);
      tot_num_records += num_records;

      const auto t1 = std::chrono::system_clock::now();
      const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
      SPDLOG_INFO(FMT_STRING("[{}:{}] processed {} records in {}s!"), chrom1.name(), chrom2.name(),
                  num_records, static_cast<double>(delta) / 1000.0);

    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format(FMT_STRING("error in while processing {}:{}: {}"),
                                           chrom1.name(), chrom2.name(), e.what()));
    } catch (...) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("An unknown error occurred while processing {}:{}"), chrom1.name(),
                      chrom2.name()));
    }
  }

  return tot_num_records;
}

static std::optional<ExpectedValues<hictk::File>> init_cis_expected_values(
    const ComputePvalConfig &c) {
  if (c.cis_only || (c.chrom1 == c.chrom2)) {
    SPDLOG_INFO(FMT_STRING("initializing expected values for cis matrices..."));
    const auto f = std::make_shared<hictk::File>(c.path_to_hic.string(), c.resolution);

    const auto bin_mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

    return {ExpectedValues<hictk::File>::cis_only(
        f,
        {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
         c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
         c.interpolation_qtile, c.interpolation_window_size},
        bin_mask)};
  }

  return {};
}

static std::size_t process_queries(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const std::optional<ExpectedValues<hictk::File>> &expected_values, const ComputePvalConfig &c) {
  if (expected_values.has_value()) {
    for (const auto &[chrom1, chrom2] : chrom_pairs) {
      try {
        if (chrom1 == chrom2) {
          std::ignore = expected_values->expected_values(chrom1);
        } else {
          std::ignore = expected_values->expected_value(chrom1, chrom2);
        }
      } catch (const std::out_of_range &) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("user-specified expected value file does not contain "
                                   "information for chromosome pair {}:{}"),
                        chrom1.name(), chrom2.name()));
      }
    }
  }

  const auto output_dir = c.output_prefix.parent_path();
  if (!output_dir.empty() && output_dir != ".") {
    std::filesystem::create_directories(output_dir);
  }

  write_chrom_sizes_to_file(hictk::File(c.path_to_hic, c.resolution).chromosomes(),
                            fmt::format(FMT_STRING("{}.chrom.sizes"), c.output_prefix), c.force);

  if (c.threads > 1) {
    BS::thread_pool tpool(conditional_static_cast<BS::concurrency_t>(c.threads));
    return process_queries_mt(tpool, chrom_pairs, expected_values, c);
  }
  return process_queries_st(chrom_pairs, expected_values, c);
}

int run_nchg_compute(const ComputePvalConfig &c) {
  const auto t0 = std::chrono::system_clock::now();

  if (c.chrom1 != "all") {
    assert(c.chrom2 != "all");
    const auto interactions_processed = run_nchg_compute_worker(c);
    const auto t1 = std::chrono::system_clock::now();
    const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    SPDLOG_INFO(FMT_STRING("[{}:{}] processed {} records in {}s!"), c.chrom1, c.chrom2,
                interactions_processed, static_cast<double>(delta) / 1000.0);
    return 0;
  }

  const hictk::File f(c.path_to_hic.string(), c.resolution);
  std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> chrom_pairs{};
  std::optional<ExpectedValues<hictk::File>> expected_values{};
  if (c.cis_only) {
    chrom_pairs = init_cis_chromosomes(f.chromosomes());
    expected_values = c.path_to_expected_values.empty()
                          ? init_cis_expected_values(c)
                          : ExpectedValues<hictk::File>::deserialize(c.path_to_expected_values);
  }
  if (c.trans_only) {
    chrom_pairs = init_trans_chromosomes(f.chromosomes());
  }

  if (c.chrom1 == "all" && !c.cis_only && !c.trans_only) {
    chrom_pairs = init_cis_chromosomes(f.chromosomes());
    const auto chrom_pairs2 = init_trans_chromosomes(f.chromosomes());
    std::copy(chrom_pairs2.begin(), chrom_pairs2.end(), std::back_inserter(chrom_pairs));
  }

  if (!expected_values.has_value() && !c.path_to_expected_values.empty()) {
    expected_values = ExpectedValues<hictk::File>::deserialize(c.path_to_expected_values);
  }

  if (expected_values.has_value() && expected_values->resolution() != f.resolution()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("mismatch in file resolution: expected values have been computed "
                               "for {}bp resolution but given Hi-C matrix has {}bp resolution"),
                    expected_values->resolution(), f.resolution()));
  }

  const auto interactions_processed = process_queries(chrom_pairs, expected_values, c);

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  SPDLOG_INFO(FMT_STRING("Processed {} records in {}s!"), interactions_processed,
              static_cast<double>(delta) / 1000.0);

  return 0;
}

}  // namespace nchg
