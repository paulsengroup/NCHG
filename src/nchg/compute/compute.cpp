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

#include <arrow/io/file.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/std.h>
#include <parquet/stream_reader.h>
#include <parquet/stream_writer.h>

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

// clang-format off
// As of HighFive 2.9.0, these headers must be included after HighFive/hictk,
// otherwise this source file fails to compile with MSVC
#include <boost/process/child.hpp>
#include <boost/process/pipe.hpp>
#include <boost/process/io.hpp>
// clang-format on

#include "nchg/common.hpp"
#include "nchg/config.hpp"
#include "nchg/io.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {

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

    std::sort(domains.begin(), domains.end());

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

  return NCHG<File>(
      f, chrom1, chrom2,
      {c.mad_max, c.min_delta, c.max_delta, c.bin_aggregation_possible_distances_cutoff,
       c.bin_aggregation_observed_distances_cutoff, c.interpolate_expected_values,
       c.interpolation_qtile, c.interpolation_window_size});
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

template <typename FilePtr>
[[nodiscard]] static std::size_t process_domains(const FilePtr &f, const ComputePvalConfig &c) {
  assert(std::filesystem::exists(c.path_to_domains));

  auto writer =
      init_parquet_file_writer<NCHGResult>(f->chromosomes(), c.output_prefix, c.force,
                                           c.compression_method, c.compression_lvl, c.threads);

  const auto domains = parse_domains(f->chromosomes(), c.path_to_domains, c.chrom1, c.chrom2);

  if (domains.empty()) {
    return 0;
  }

  const auto nchg = init_nchg(f, c);

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

template <typename FilePtr>
[[nodiscard]] static std::size_t process_one_chromosome_pair(const FilePtr &f,
                                                             const ComputePvalConfig &c) {
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);
  auto nchg = init_nchg(f, c);

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
      "--verbosity",
      "2",
      "--min-delta",
      fmt::to_string(c.min_delta),
      "--max-delta",
      fmt::to_string(c.max_delta),
      "--resolution",
      fmt::to_string(c.resolution),
      c.path_to_hic.string(),
      output_path.string(),
  };

  if (!c.path_to_domains.empty()) {
    args.emplace_back("--domains");
    args.emplace_back(c.path_to_domains.string());
  }

  if (!c.path_to_expected_values.empty()) {
    args.emplace_back("--expected-values");
    args.emplace_back(c.path_to_expected_values.string());
  }

  if (c.force) {
    args.emplace_back("--force");
  }

  boost::process::child proc(
      c.exec.string(), args,
      boost::process::std_in<boost::process::null, boost::process::std_out> boost::process::null);
  SPDLOG_DEBUG(FMT_STRING("spawned worker process {}..."), proc.id());
  if (!proc.running()) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to spawn worker process: {} {}"),
                                         c.exec.string(), fmt::join(args, " ")));
  }

  return proc;
}

static std::size_t process_queries_mt(
    BS::thread_pool &tpool,
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const ComputePvalConfig &c) {
  std::atomic<bool> early_return{false};

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
              fmt::format(FMT_STRING("{}.{}.{}.parquet"), c.output_prefix.string(), chrom1.name(),
                          chrom2.name());

          auto proc = spawn_compute_process(c, output_path, chrom1, chrom2);
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

static std::size_t process_queries_st(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const ComputePvalConfig &c) {
  std::size_t tot_num_records = 0;
  for (const auto &pair : chrom_pairs) {
    const auto &chrom1 = pair.first;
    const auto &chrom2 = pair.second;

    try {
      auto config = c;
      config.chrom1 = chrom1.name();
      config.chrom2 = chrom2.name();

      config.output_prefix = fmt::format(FMT_STRING("{}.{}.{}.parquet"), c.output_prefix.string(),
                                         chrom1.name(), chrom2.name());

      const auto t0 = std::chrono::system_clock::now();
      const auto num_records = run_nchg_compute_worker(config);
      tot_num_records += num_records;

      const auto t1 = std::chrono::system_clock::now();
      const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
      SPDLOG_INFO(FMT_STRING("Processed {} records in {}s!"), num_records,
                  static_cast<double>(delta) / 1000.0);

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

static std::size_t process_queries(
    const std::vector<std::pair<hictk::Chromosome, hictk::Chromosome>> &chrom_pairs,
    const ComputePvalConfig &c) {
  const auto output_dir = c.output_prefix.parent_path();
  if (!output_dir.empty() && output_dir != ".") {
    std::filesystem::create_directories(output_dir);
  }

  write_chrom_sizes_to_file(hictk::File(c.path_to_hic, c.resolution).chromosomes(),
                            fmt::format(FMT_STRING("{}.chrom.sizes"), c.output_prefix), c.force);

  if (c.threads > 1) {
    BS::thread_pool tpool(conditional_static_cast<BS::concurrency_t>(c.threads));
    return process_queries_mt(tpool, chrom_pairs, c);
  }
  return process_queries_st(chrom_pairs, c);
}

int run_nchg_compute(const ComputePvalConfig &c) {
  const auto t0 = std::chrono::system_clock::now();

  if (c.chrom1 != "all") {
    assert(c.chrom2 != "all");
    const auto interactions_processed = run_nchg_compute_worker(c);
    const auto t1 = std::chrono::system_clock::now();
    const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    SPDLOG_INFO(FMT_STRING("Processed {} records in {}s!"), interactions_processed,
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

  const auto interactions_processed = process_queries(chrom_pairs, c);

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  SPDLOG_INFO(FMT_STRING("Processed {} records in {}s!"), interactions_processed,
              static_cast<double>(delta) / 1000.0);

  return 0;
}

}  // namespace nchg
