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

#include "nchg/file_store.hpp"

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
#include "nchg/parquet_stats_file_writer.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

using ChromPair = std::pair<hictk::Chromosome, hictk::Chromosome>;
using BG2Domain = std::pair<hictk::GenomicInterval, hictk::GenomicInterval>;

using BG2DomainSet = phmap::flat_hash_set<BG2Domain>;
using ChromosomePairs = phmap::btree_set<ChromPair>;

class BG2Domains {
  std::vector<BG2Domain> _domains;
  std::unique_ptr<std::mutex> _mtx{};

  [[nodiscard]] static std::vector<BG2Domain> parse_domains(
      const hictk::Reference &chroms, const std::filesystem::path &path, bool keep_cis,
      bool keep_trans, const std::optional<hictk::Chromosome> &chrom1,
      const std::optional<hictk::Chromosome> &chrom2) {
    SPDLOG_INFO("reading domains from file \"{}\"...", path.string());
    const auto t0 = std::chrono::steady_clock::now();

    std::ifstream ifs{};
    ifs.exceptions(ifs.exceptions() | std::ios::badbit | std::ios::failbit);

    std::size_t i = 1;
    std::string buffer;
    std::size_t domains_dropped{};
    std::size_t duplicate_domains{};
    std::size_t dropped_cis{};
    std::size_t dropped_trans{};
    BG2DomainSet domain_set{};

    if (chrom1.has_value()) {
      assert(chrom2.has_value());
      assert(chroms.contains(*chrom1));
      assert(chroms.contains(*chrom2));
    }

    try {
      ifs.open(path);

      for (; std::getline(ifs, buffer); ++i) {
        if (buffer.empty()) {
          continue;
        }

        if (buffer.back() == '\r') {
          buffer.resize(buffer.size() - 1);
        }

        const auto record = truncate_record<6>(buffer);
        const auto domain1 = truncate_record<3>(record);
        const auto domain2 = truncate_record<3>(record.substr(domain1.size() + 1));

        const auto chrom1_parsed = truncate_record<1>(domain1);
        const auto chrom2_parsed = truncate_record<1>(domain2);

        if (!chrom1.has_value()) {
          assert(!chrom2.has_value());
          if (!chroms.contains(chrom1_parsed) || !chroms.contains(chrom2_parsed)) {
            ++domains_dropped;
            continue;
          }

          if (!keep_cis && chrom1_parsed == chrom2_parsed) {
            ++dropped_cis;
            continue;
          }

          if (!keep_trans && chrom1_parsed != chrom2_parsed) {
            ++dropped_trans;
            continue;
          }
        }

        if ((chrom1.has_value() && *chrom1 != chrom1_parsed) ||
            (chrom2.has_value() && *chrom2 != chrom2_parsed)) {
          ++domains_dropped;
          continue;
        }

        auto domain = std::make_pair(hictk::GenomicInterval::parse_bed(chroms, domain1),
                                     hictk::GenomicInterval::parse_bed(chroms, domain2));

        if (domain.first > domain.second) {
          throw std::runtime_error(
              fmt::format("domains cannot overlap with the lower triangular matrix: offending "
                          "domain {:ucsc}; {:ucsc}",
                          domain.first, domain.second));
        }

        const auto &[_, inserted] = domain_set.emplace(std::move(domain));
        duplicate_domains += static_cast<std::size_t>(!inserted);
      }
    } catch (const std::exception &e) {
      if (!ifs.eof()) {
        throw std::runtime_error(
            fmt::format("found an invalid record at line {} of file {}: {}", i, path, e.what()));
      }
    } catch (...) {
      throw std::runtime_error(
          fmt::format("found an invalid record at line {} of file {}: unknown error", i, path));
    }

    if (duplicate_domains != 0) {
      SPDLOG_WARN("found {} duplicate domain(s)", domains_dropped);
    }
    if (domains_dropped != 0) {
      if (chrom1.has_value()) {
        assert(chrom2.has_value());
        SPDLOG_WARN(
            "[{}:{}]: {} domain(s) were dropped because they did not map to the selected "
            "chromosomes",
            chrom1->name(), chrom2->name(), domains_dropped);
      } else {
        SPDLOG_WARN("{} domain(s) were dropped because they did not map to any known chromosome",
                    domains_dropped);
      }
    }
    if (dropped_cis != 0) {
      SPDLOG_WARN(
          "{} domain(s) were dropped because they overlapped with the cis area of the interaction "
          "map",
          dropped_cis);
    }
    if (dropped_trans != 0) {
      SPDLOG_WARN(
          "{} domain(s) were dropped because they overlapped with the trans area of the "
          "interaction map",
          dropped_trans);
    }

    if (domain_set.empty()) {
      throw std::runtime_error(
          fmt::format("unable to parse any domain from file \"{}\"", path.string()));
    }

    std::vector<BG2Domain> domains(domain_set.size());
    std::ranges::move(domain_set, domains.begin());
    {
      BG2DomainSet tmp{};
      std::swap(domain_set, tmp);
    }

    // We want to domains to be sorted in a tiled fashion:
    // e.g. all chr1:chr1 domains should precede all chr1:chr2 domains.
    //      Within tiles, domains should be sorted by their genomic coordinates
    std::ranges::sort(domains, [](const auto &domain1, const auto &domain2) {
      const auto &c1 = domain1.first.chrom();
      const auto &c2 = domain2.first.chrom();
      const auto &c3 = domain1.second.chrom();
      const auto &c4 = domain2.second.chrom();

      if (c1 != c2) {
        return c1 < c2;
      }
      if (c3 != c4) {
        return c3 < c4;
      }

      auto pos1 = domain1.first.start();
      auto pos2 = domain2.first.start();

      if (pos1 != pos2) {
        return pos1 < pos2;
      }

      auto pos3 = domain1.second.start();
      auto pos4 = domain2.second.start();

      if (pos3 != pos4) {
        return pos3 < pos4;
      }

      pos1 = domain1.first.end();
      pos2 = domain2.first.end();

      if (pos1 != pos2) {
        return pos1 < pos2;
      }

      pos3 = domain1.second.end();
      pos4 = domain2.second.end();

      return pos3 < pos4;
    });

    const auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO("read {} domains from \"{}\" in {}", domains.size(), path.string(),
                format_duration(t1 - t0));
    return domains;
  }
  [[nodiscard]] auto select_domains(const hictk::Chromosome &chrom1,
                                    const hictk::Chromosome &chrom2) noexcept {
    SPDLOG_DEBUG("[{}:{}]: selecting domains...", chrom1.name(), chrom2.name());

    // The positions here do not matter: see implementation of the comparison operator below
    const auto query =
        std::make_pair(hictk::GenomicInterval{chrom1, 0, 1}, hictk::GenomicInterval{chrom2, 0, 1});
    return std::equal_range(_domains.begin(), _domains.end(), query,
                            [&](const auto &domain1, const auto &domain2) {
                              const auto &c1 = domain1.first.chrom();
                              const auto &c2 = domain2.first.chrom();
                              const auto &c3 = domain1.second.chrom();
                              const auto &c4 = domain2.second.chrom();

                              if (c1 == c2) {
                                return c3 < c4;
                              }
                              return c1 < c2;
                            });
  }

 public:
  using iterator = std::vector<BG2Domain>::iterator;
  using const_iterator = std::vector<BG2Domain>::const_iterator;

  BG2Domains() : _mtx(std::make_unique<std::mutex>()) {}
  BG2Domains(const hictk::Reference &chroms, const std::filesystem::path &path, bool keep_cis,
             bool keep_trans, const std::optional<hictk::Chromosome> &chrom1 = {},
             const std::optional<hictk::Chromosome> &chrom2 = {})
      : _domains(parse_domains(chroms, path, keep_cis, keep_trans, chrom1, chrom2)),
        _mtx(std::make_unique<std::mutex>()) {}

  [[nodiscard]] iterator begin() noexcept { return _domains.begin(); }
  [[nodiscard]] iterator end() noexcept { return _domains.end(); }

  [[nodiscard]] const_iterator begin() const noexcept { return _domains.begin(); }
  [[nodiscard]] const_iterator end() const noexcept { return _domains.end(); }

  [[nodiscard]] const_iterator cbegin() const noexcept { return _domains.cbegin(); }
  [[nodiscard]] const_iterator cend() const noexcept { return _domains.cend(); }

  [[nodiscard]] std::vector<BG2Domain> extract(const hictk::Chromosome &chrom1,
                                               const hictk::Chromosome &chrom2) {
    assert(_mtx);
    SPDLOG_DEBUG("[{}:{}] extracting domains...", chrom1.name(), chrom2.name());

    [[maybe_unused]] const auto lck = std::scoped_lock(*_mtx);
    auto [first, last] = select_domains(chrom1, chrom2);
    std::vector domains(first, last);
    _domains.erase(first, last);
    _domains.shrink_to_fit();
    SPDLOG_DEBUG("[{}:{}] extracted {} domains!", chrom1.name(), chrom2.name(), domains.size());
    return domains;
  }

  [[nodiscard]] std::filesystem::path extact_and_write_to_file(
      const std::filesystem::path &dest_dir, const hictk::Chromosome &chrom1,
      const hictk::Chromosome &chrom2, bool force) {
    const auto dest = dest_dir / fmt::format("domains.{}.{}.bedpe", chrom1.name(), chrom2.name());
    SPDLOG_DEBUG("[{}:{}] writing domains to temporary file \"{}\"...", chrom1.name(),
                 chrom2.name(), dest.string());

    const auto t0 = std::chrono::steady_clock::now();

    if (!force && std::filesystem::exists(dest)) {
      throw std::runtime_error(
          fmt::format("refusing to overwrite existing temporary file: \"{}\". "
                      "Pass --force to overwrite.",
                      dest));
    }

    std::filesystem::remove(dest);  // NOLINT;

    [[maybe_unused]] std::size_t domains_processed{};
    std::ofstream fs{};
    fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

    try {
      if (dest_dir.empty() && !std::filesystem::exists(dest)) {
        std::filesystem::create_directories(dest_dir);  // NOLINT
      }
#ifdef __cpp_lib_ios_noreplace
      fs.open(dest, std::ios::out | std::ios::trunc | std::ios::noreplace);
#else
      fs.open(dest, std::ios::out | std::ios::trunc);
#endif

      const auto selected_domains = extract(chrom1, chrom2);
      if (selected_domains.empty()) {
        SPDLOG_WARN("[{}:{}]: no domains were selected!", chrom1.name(), chrom2.name());
      } else if (selected_domains.size() == 1) {
        SPDLOG_DEBUG("[{}:{}] selected a single domain ({:ucsc}; {:ucsc})...", chrom1.name(),
                     chrom2.name(), selected_domains.front().first,
                     selected_domains.front().second);
      } else {
        SPDLOG_DEBUG("[{}:{}] selected {} domains ({:ucsc}; {:ucsc} ... {:ucsc}; {:ucsc})...",
                     chrom1.name(), chrom2.name(), selected_domains.size(),
                     selected_domains.front().first, selected_domains.front().second,
                     selected_domains.back().first, selected_domains.back().second);
      }
      for (const auto &domain : selected_domains) {
        fmt::print(fs, "{:bed}\t{:bed}\n", domain.first, domain.second);
        ++domains_processed;
      }

    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format("failed to write domains for {}:{} to temporary file \"{}\": {}",
                      chrom1.name(), chrom2.name(), dest.string(), e.what()));
    } catch (...) {
      throw std::runtime_error(
          fmt::format("failed to write domains for {}:{} to temporary file \"{}\": unknown error",
                      chrom1.name(), chrom2.name(), dest.string()));
    }

    const auto t1 = std::chrono::steady_clock::now();

    SPDLOG_DEBUG("[{}:{}] written {} domains to \"{}\" in {}", chrom1.name(), chrom2.name(),
                 domains_processed, dest.string(), format_duration(t1 - t0));
    return dest;
  }
};

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
#ifdef __cpp_lib_ios_noreplace
    fs.open(path, std::ios::out | std::ios::noreplace);
#else
    fs.open(path);
#endif

    for (const auto &chrom : chroms) {
      fmt::print(fs, "{}\t{}\n", chrom.name(), chrom.size());
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to write chromosomes to file {}: {}", path, e.what()));
  }
}

[[nodiscard]] static std::size_t process_domains(
    const std::shared_ptr<const hictk::File> &f, BG2Domains &domains,
    const std::optional<ExpectedValues> &expected_values, const ComputePvalConfig &c) {
  assert(std::filesystem::exists(c.path_to_domains));
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(!c.output_path.empty());

  SPDLOG_INFO("[{}:{}] begin processing domains from {}...", *c.chrom1, *c.chrom2,
              c.path_to_domains);

  ParquetStatsFileWriter writer(f->chromosomes(), c.output_path, c.force, c.compression_method,
                                c.compression_lvl, c.threads);

  const auto selected_domains =
      domains.extract(f->chromosomes().at(*c.chrom1), f->chromosomes().at(*c.chrom2));
  if (selected_domains.empty()) {
    writer.finalize<NCHGResult>();
    return 0;
  }

  const auto nchg = init_nchg(f, expected_values, c);

  std::size_t num_records = 0;
  for (const auto &domain : selected_domains) {
    const auto s = nchg.compute(domain.first, domain.second, c.bad_bin_fraction);

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

  SPDLOG_INFO("[{}:{}] begin processing interactions...", *c.chrom1, *c.chrom2);

  const auto &chrom1 = f->chromosomes().at(*c.chrom1);
  const auto &chrom2 = f->chromosomes().at(*c.chrom2);
  const auto nchg = init_nchg(f, expected_values, c);

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

[[nodiscard]] static std::size_t run_nchg_compute_worker(
    const ComputePvalConfig &c, std::optional<BG2Domains> &domains,
    const std::optional<ExpectedValues> &expected_values = {}) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());

  const auto f = std::make_shared<const hictk::File>(c.path_to_hic.string(), c.resolution);

  if (domains.has_value()) {
    return process_domains(f, *domains, expected_values, c);
  }

  if (!c.path_to_domains.empty()) {
    const auto chrom1 = std::make_optional(f->chromosomes().at(*c.chrom1));
    const auto chrom2 = std::make_optional(f->chromosomes().at(*c.chrom2));

    BG2Domains domains_(f->chromosomes(), c.path_to_domains, true, true, chrom1, chrom2);
    return process_domains(f, domains_, expected_values, c);
  }

  return process_one_chromosome_pair(f, expected_values, c);
}

[[nodiscard]] static ChromosomePairs init_cis_chromosomes(
    const hictk::Reference &chroms, const std::optional<BG2Domains> &domains) {
  ChromosomePairs buffer{};

  if (domains.has_value()) {
    for (const auto &[domain1, domain2] : *domains) {
      assert(chroms.contains(domain1.chrom().name()));
      assert(chroms.contains(domain2.chrom().name()));

      if (domain1.chrom() == domain2.chrom()) {
        buffer.emplace(domain1.chrom(), domain2.chrom());
      }
    }
  } else {
    for (const auto &chrom : chroms) {
      if (chrom.is_all()) [[unlikely]] {
        continue;
      }
      buffer.emplace(chrom, chrom);
    }
  }

  return buffer;
}

[[nodiscard]] static ChromosomePairs init_trans_chromosomes(
    const hictk::Reference &chroms, const std::optional<BG2Domains> &domains) {
  ChromosomePairs buffer{};

  if (domains.has_value()) {
    for (const auto &[domain1, domain2] : *domains) {
      assert(chroms.contains(domain1.chrom().name()));
      assert(chroms.contains(domain2.chrom().name()));

      if (domain1.chrom() != domain2.chrom()) {
        buffer.emplace(domain1.chrom(), domain2.chrom());
      }
    }
  } else {
    for (std::uint32_t chrom1_id = 0; chrom1_id < chroms.size(); ++chrom1_id) {
      const auto &chrom1 = chroms.at(chrom1_id);
      if (chrom1.is_all()) [[unlikely]] {
        continue;
      }
      for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < chroms.size(); ++chrom2_id) {
        const auto &chrom2 = chroms.at(chrom2_id);
        if (chrom2.is_all()) [[unlikely]] {
          continue;
        }
        buffer.emplace(chrom1, chrom2);
      }
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

[[nodiscard]] static std::filesystem::path generate_report_file_name(
    const std::filesystem::path &output_prefix) {
  assert(!output_prefix.empty());
  return fmt::format("{}.json", output_prefix);
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

[[nodiscard]] static std::size_t worker_fx(FileStore &file_store, const hictk::Chromosome &chrom1,
                                           const hictk::Chromosome &chrom2,
                                           std::optional<BG2Domains> &domains,
                                           const ComputePvalConfig &config,
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
      domain_file = domains->extact_and_write_to_file(child_config.tmpdir, chrom1, chrom2,
                                                      child_config.force);
      child_config.path_to_domains = domain_file;
    } else {
      assert(child_config.path_to_domains.empty());
    }

    if (!trans_expected_values_avail && chrom1 != chrom2) {
      child_config.path_to_expected_values.clear();
    }

    auto proc = spawn_compute_process(child_config, chrom1, chrom2);
    proc.wait();

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

static std::size_t process_queries_mt(BS::thread_pool &tpool, FileStore &file_store,
                                      const ChromosomePairs &chrom_pairs,
                                      std::optional<BG2Domains> &domains,
                                      const std::optional<ExpectedValues> &expected_values,
                                      const ComputePvalConfig &c) {
  std::atomic<bool> early_return{false};

  const auto user_provided_expected_values = !c.path_to_expected_values.empty();
  const auto base_config = init_base_config(c, expected_values);

  BS::multi_future<std::size_t> workers(chrom_pairs.size());
  auto it = chrom_pairs.begin();
  for (std::size_t i = 0; i < workers.size(); ++i) {
    workers[i] = tpool.submit_task([&, i, chrom_pair = *it++] {
      const auto &[chrom1, chrom2] = chrom_pair;
      SPDLOG_DEBUG("submitting task {}/{} ({}:{})...", i + 1, workers.size(), chrom1.name(),
                   chrom2.name());
      return worker_fx(file_store, chrom1, chrom2, domains, base_config,
                       user_provided_expected_values, early_return);
    });
  }

  std::size_t num_records = 0;
  for (const auto result : workers.get()) {
    num_records += result;
  }

  return num_records;
}

static std::size_t process_queries_st(FileStore &file_store, const ChromosomePairs &chrom_pairs,
                                      std::optional<BG2Domains> &domains,
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
                                   ? run_nchg_compute_worker(config, domains, expected_values)
                                   : run_nchg_compute_worker(config, domains);
      tot_num_records += num_records;

      file_store.register_file(c.output_path);

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

static std::size_t process_queries(FileStore &file_store, const ChromosomePairs &chrom_pairs,
                                   std::optional<BG2Domains> &domains,
                                   const std::optional<ExpectedValues> &expected_values,
                                   const ComputePvalConfig &c) {
  const std::filesystem::path chrom_sizes_path{generate_chrom_sizes_file_name(c.output_prefix)};
  write_chrom_sizes_to_file(hictk::File(c.path_to_hic, c.resolution).chromosomes(),
                            chrom_sizes_path, c.force);
  file_store.register_file(chrom_sizes_path);

  if (c.threads > 1) {
    auto num_threads = conditional_static_cast<BS::concurrency_t>(c.threads);
    if (c.threads > chrom_pairs.size()) {
      num_threads = conditional_static_cast<BS::concurrency_t>(chrom_pairs.size());
      SPDLOG_WARN(
          "number of threads specified through --threads exceeds the number of chromosome pairs "
          "to "
          "be processed: limiting concurrency to {} thread(s)",
          num_threads);
    }
    BS::thread_pool tpool(num_threads);
    return process_queries_mt(tpool, file_store, chrom_pairs, domains, expected_values, c);
  }
  return process_queries_st(file_store, chrom_pairs, domains, expected_values, c);
}

static void process_file_collisions(const std::filesystem::path &output_prefix,
                                    const ChromosomePairs &chrom_pairs, bool force) {
  std::vector<std::string> collisions{};
  std::size_t num_collisions = 0;

  constexpr std::size_t max_collisions_reported = 10;

  if (const auto report_file = generate_report_file_name(output_prefix); force) {
    const auto removed = std::filesystem::remove(report_file);  // NOLINT
    if (removed) {
      SPDLOG_DEBUG("file \"{}\" has been deleted", report_file.string());
    }
  } else if (std::filesystem::exists(report_file)) {
    collisions.emplace_back(report_file.string());
    ++num_collisions;
  }

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

static void validate_expected_values(const ExpectedValues &expected_values,
                                     const std::filesystem::path &path_to_expected_values,
                                     const ChromosomePairs &chrom_pairs, std::uint32_t resolution) {
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
    std::unique_ptr<FileStore> file_store;
    ChromosomePairs chrom_pairs{};
    std::optional<ExpectedValues> expected_values{};
    std::optional<BG2Domains> domains{};
  };

  Plan plan{};

  const hictk::File f(c.path_to_hic.string(), c.resolution);

  if (f.bins().type() != hictk::BinTable::Type::fixed) {
    throw std::runtime_error("only file with uniform bin sizes are supported.");
  }

  if (!c.path_to_domains.empty()) {
    plan.domains.emplace(f.chromosomes(), c.path_to_domains, c.compute_cis, c.compute_trans);
  }

  if (c.compute_cis) {
    plan.chrom_pairs = init_cis_chromosomes(f.chromosomes(), plan.domains);
  }

  if (c.compute_trans) {
    std::ranges::move(init_trans_chromosomes(f.chromosomes(), plan.domains),
                      std::inserter(plan.chrom_pairs, plan.chrom_pairs.end()));
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

  const auto root_dir = c.output_prefix.parent_path();

  if (!root_dir.empty() && !std::filesystem::exists(root_dir)) {
    std::filesystem::create_directories(root_dir);  // NOLINT
  }
  plan.file_store = std::make_unique<FileStore>(
      root_dir, false, fmt::format("{}.json", c.output_prefix.filename().string()));

  return plan;
}

int run_command(const ComputePvalConfig &c) {
  try {
    const auto t0 = std::chrono::steady_clock::now();

    if (c.chrom1.has_value()) {
      assert(c.chrom2.has_value());
      assert(!c.output_path.empty());
      std::optional<BG2Domains> placeholder{};
      const auto interactions_processed = run_nchg_compute_worker(c, placeholder);
      const auto t1 = std::chrono::steady_clock::now();
      SPDLOG_INFO("[{}:{}] processed {} records in {}!", *c.chrom1, *c.chrom2,
                  interactions_processed, format_duration(t1 - t0));

      SPDLOG_INFO("records for {}:{} have been written to file \"{}\"", *c.chrom1, *c.chrom2,
                  c.output_path.string());
      return 0;
    }

    auto [file_store, chrom_pairs, expected_values, domains] = generate_execution_plan(c);

    const auto interactions_processed =
        process_queries(*file_store, chrom_pairs, domains, expected_values, c);
    file_store->finalize();
    std::filesystem::remove_all(c.tmpdir);  // NOLINT

    const auto t1 = std::chrono::steady_clock::now();
    if (interactions_processed == 0) {
      SPDLOG_WARN("no records have been processed. Is this intended?");
    } else {
      SPDLOG_INFO("processed {} records in {}!", interactions_processed, format_duration(t1 - t0));
    }

    SPDLOG_INFO("{} new file(s) have been created under prefix \"{}\"", chrom_pairs.size() + 2,
                c.output_prefix);

    return 0;
  } catch (...) {
    std::error_code ec{};
    std::filesystem::remove_all(c.tmpdir, ec);  // NOLINT

    if (!static_cast<bool>(ec) && std::filesystem::exists(c.tmpdir)) {
      SPDLOG_WARN("failed to remove temporary folder \"{}\"! Please remove the folder manually",
                  c.tmpdir.string());
    }
    throw;
  }
}

}  // namespace nchg
