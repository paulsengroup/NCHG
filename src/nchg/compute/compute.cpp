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
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
#include <parquet/stream_reader.h>
#include <parquet/stream_writer.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <glaze/glaze.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "nchg/tools/tmpdir.hpp"

// clang-format off
// As of HighFive 2.9.0, these headers must be included after HighFive/hictk,
// otherwise this source file fails to compile with MSVC
#include <boost/asio/io_context.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/process/v2/environment.hpp>
#include <boost/process/v2/process.hpp>
#include <boost/process/v2/stdio.hpp>
// clang-format on

#include "nchg/common.hpp"
#include "nchg/file_store.hpp"
#include "nchg/genomic_domains.hpp"
#include "nchg/nchg.hpp"
#include "nchg/parquet_stats_file_writer.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace glz {

template <>
struct from<BEVE, spdlog::log_clock::time_point> {
  template <auto Opts>
  static void op(spdlog::log_clock::time_point &timestamp, auto &&...args) {
    using BuffT = decltype(timestamp.time_since_epoch().count());

    BuffT buff{};
    parse<BEVE>::op<Opts>(buff, args...);
    timestamp = spdlog::log_clock::time_point(spdlog::log_clock::duration{buff});
  }
};

template <>
struct from<BEVE, spdlog::string_view_t> {
  template <auto Opts>
  static void op(spdlog::string_view_t &payload, auto &&...args) {
    std::string_view buff{};
    parse<BEVE>::op<Opts>(buff, args...);
    payload = spdlog::string_view_t{buff.data(), buff.size()};
  }
};

template <>
struct to<BEVE, spdlog::log_clock::time_point> {
  template <auto Opts>
  static void op(const spdlog::log_clock::time_point &timestamp, auto &&...args) {
    serialize<BEVE>::op<Opts>(timestamp.time_since_epoch().count(), args...);
  }
};

template <>
struct to<BEVE, spdlog::string_view_t> {
  template <auto Opts>
  static void op(spdlog::string_view_t s, auto &&...args) {
    serialize<BEVE>::op<Opts>(std::string_view{s.data(), s.size()}, args...);
  }
};

template <>
struct meta<spdlog::details::log_msg> {
  using T = spdlog::details::log_msg;
  static constexpr auto value = object("t", &T::time, "l", &T::level, "p", &T::payload);
};

}  // namespace glz

namespace nchg {

using DomainAggregationStrategy = ComputePvalConfig::DomainAggregationStrategy;
using ChromPair = std::pair<hictk::Chromosome, hictk::Chromosome>;
using ChromosomePairs = phmap::btree_set<ChromPair>;

struct GIHasher {
  using is_transparent = void;
  static std::size_t operator()(const std::shared_ptr<const hictk::GenomicInterval> &gi) {
    assert(!!gi);
    return operator()(*gi);
  }

  static std::size_t operator()(const hictk::GenomicInterval &gi) {
    return std::hash<hictk::GenomicInterval>{}(gi);
  }
};

struct GIEqOperator {
  using is_transparent = void;
  static bool operator()(const std::shared_ptr<const hictk::GenomicInterval> &gi1,
                         const std::shared_ptr<const hictk::GenomicInterval> &gi2) {
    assert(!!gi1);
    assert(!!gi2);

    return *gi1 == *gi2;
  }

  static bool operator()(const hictk::GenomicInterval &gi1,
                         const std::shared_ptr<const hictk::GenomicInterval> &gi2) {
    assert(!!gi2);

    return gi1 == *gi2;
  }

  static bool operator()(const std::shared_ptr<const hictk::GenomicInterval> &gi1,
                         const hictk::GenomicInterval &gi2) {
    assert(!!gi1);

    return *gi1 == gi2;
  }
};

[[nodiscard]] static GenomicDomains parse_domains(
    const hictk::Reference &chroms, const std::filesystem::path &path, bool keep_cis,
    bool keep_trans, const std::optional<hictk::Chromosome> &chrom1 = {},
    const std::optional<hictk::Chromosome> &chrom2 = {}) {
  SPDLOG_INFO("reading domains from file \"{}\"...", path.string());
  const auto t0 = std::chrono::steady_clock::now();

  std::ifstream ifs{};
  ifs.exceptions(ifs.exceptions() | std::ios::badbit | std::ios::failbit);

  std::size_t i = 1;
  std::string buffer;
  std::size_t domains_parsed{};
  std::size_t domains_dropped{};
  std::size_t duplicate_domains{};
  std::size_t dropped_cis{};
  std::size_t dropped_trans{};

  phmap::flat_hash_set<std::shared_ptr<const hictk::GenomicInterval>, GIHasher, GIEqOperator>
      intervals{};
  phmap::flat_hash_set<BEDPE> domains{};

  if (chrom1.has_value()) {
    assert(chrom2.has_value());
    assert(chroms.contains(*chrom1));
    assert(chroms.contains(*chrom2));
  }

  auto fetch_or_insert_gi = [&](std::string_view domain) {
    auto gi = hictk::GenomicInterval::parse_bed(chroms, domain);
    const auto match = intervals.find(gi);
    if (match != intervals.end()) {
      return *match;
    }

    const auto [it, _] =
        intervals.emplace(std::make_shared<const hictk::GenomicInterval>(std::move(gi)));
    return *it;
  };

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

      ++domains_parsed;

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

      auto gi1 = std::make_shared<const hictk::GenomicInterval>(
          hictk::GenomicInterval::parse_bed(chroms, domain1));
      auto gi2 = domain1 == domain2 ? gi1
                                    : std::make_shared<const hictk::GenomicInterval>(
                                          hictk::GenomicInterval::parse_bed(chroms, domain2));

      BEDPE domain{fetch_or_insert_gi(domain1), fetch_or_insert_gi(domain2)};

      if (domain.range1() > domain.range2()) {
        throw std::runtime_error(
            fmt::format("domains cannot overlap with the lower triangular matrix: offending "
                        "domain {:ucsc}; {:ucsc}",
                        domain.range1(), domain.range2()));
      }

      const auto &[_, inserted] = domains.emplace(std::move(domain));
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
    SPDLOG_WARN("found {} duplicate domain(s)", duplicate_domains);
  }
  if (domains_dropped != 0) {
    if (chrom1.has_value()) {
      assert(chrom2.has_value());
      SPDLOG_DEBUG(
          "[{}:{}]: {}/{} domain(s) were dropped because they did not map to the specified "
          "chromosome(s)",
          chrom1->name(), chrom2->name(), domains_dropped, domains_parsed);
    } else {
      SPDLOG_WARN("{}/{} domain(s) were dropped because they did not map to any known chromosome",
                  domains_dropped, domains_parsed);
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

  if (domains.empty()) {
    throw std::runtime_error(
        fmt::format("unable to parse any domain from file \"{}\"", path.string()));
  }

  GenomicDomains domains_{std::vector<BEDPE>{std::make_move_iterator(domains.begin()),
                                             std::make_move_iterator(domains.end())},
                          true};

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("read {} domains from \"{}\" in {}", domains_.size(), path.string(),
              format_duration(t1 - t0));
  return domains_;
}

[[nodiscard]] static GenomicDomains parse_domains(const hictk::Reference &chroms,
                                                  const std::filesystem::path &path, bool keep_cis,
                                                  bool keep_trans,
                                                  const std::optional<std::string> &chrom1 = {},
                                                  const std::optional<std::string> &chrom2 = {}) {
  if (!chrom1.has_value()) {
    assert(!chrom2.has_value());
    return parse_domains(chroms, path, keep_cis, keep_trans, std::optional<hictk::Chromosome>{},
                         std::optional<hictk::Chromosome>{});
  }

  assert(chrom2.has_value());
  assert(chroms.contains(*chrom1));
  assert(chroms.contains(*chrom2));

  return parse_domains(chroms, path, keep_cis, keep_trans, chroms.at(*chrom1), chroms.at(*chrom2));
}

[[nodiscard]] static std::filesystem::path write_domains_to_file(
    const GenomicDomains &domains, const std::filesystem::path &dest_dir,
    const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2, bool force) {
  const auto dest = dest_dir / fmt::format("domains.{}.{}.bedpe", chrom1.name(), chrom2.name());
  SPDLOG_DEBUG("[{}:{}]: writing domains to temporary file \"{}\"...", chrom1.name(), chrom2.name(),
               dest.string());

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

    const auto selected_domains = domains.fetch<std::uint8_t>(chrom1, chrom2).to_vector();
    if (selected_domains.empty()) {
      SPDLOG_WARN("[{}:{}]: no domains were selected!", chrom1.name(), chrom2.name());
    } else if (selected_domains.size() == 1) {
      SPDLOG_DEBUG("[{}:{}]: selected a single domain ({:ucsc}; {:ucsc})...", chrom1.name(),
                   chrom2.name(), selected_domains.front().first.range1(),
                   selected_domains.front().first.range2());
    } else {
      SPDLOG_DEBUG("[{}:{}]: selected {} domains ({:ucsc}; {:ucsc} ... {:ucsc}; {:ucsc})...",
                   chrom1.name(), chrom2.name(), selected_domains.size(),
                   selected_domains.front().first.range1(), selected_domains.front().first.range2(),
                   selected_domains.back().first.range1(), selected_domains.back().first.range2());
    }

    for (const auto &[domain, _] : selected_domains) {
      fmt::print(fs, "{:bed}\t{:bed}\n", domain.range1(), domain.range2());
      ++domains_processed;
    }

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("failed to write domains for {}:{} to temporary file \"{}\": {}", chrom1.name(),
                    chrom2.name(), dest.string(), e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to write domains for {}:{} to temporary file \"{}\": unknown error",
                    chrom1.name(), chrom2.name(), dest.string()));
  }

  const auto t1 = std::chrono::steady_clock::now();

  SPDLOG_DEBUG("[{}:{}]: written {} domains to \"{}\" in {}", chrom1.name(), chrom2.name(),
               domains_processed, dest.string(), format_duration(t1 - t0));
  return dest;
}

class ProcessContext {
  using LockedContext = std::pair<std::unique_lock<std::mutex>, boost::asio::io_context *>;
  boost::asio::io_context _ctx;
  std::mutex _mtx;

 public:
  ProcessContext() = default;

  [[nodiscard]] LockedContext operator()() { return std::make_pair(std::unique_lock{_mtx}, &_ctx); }
};

[[nodiscard]] static auto init_nchg(const std::shared_ptr<const hictk::File> &f,
                                    const std::optional<ExpectedValues> &expected_values,
                                    const ComputePvalConfig &c) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());

  struct Result {
    NCHG nchg;
    std::shared_ptr<const std::vector<bool>> bin1_mask{};
    std::shared_ptr<const std::vector<bool>> bin2_mask{};
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

[[nodiscard]] static std::uint64_t map_interactions_to_domains_one_pass(
    const hictk::File &f, GenomicDomainsIndexed<std::uint64_t> &obs_domains,
    GenomicDomainsIndexed<double> &exp_domains, const ExpectedMatrixStats &expected_matrix,
    std::uint64_t min_delta, std::uint64_t max_delta, const std::vector<bool> &bin1_mask,
    const std::vector<bool> &bin2_mask) {
  const auto &chrom1 = obs_domains.chrom1().name();
  const auto &chrom2 = obs_domains.chrom2().name();

  SPDLOG_DEBUG("[{}:{}]: mapping interactions to {} genomic domains using the one-pass strategy...",
               chrom1, chrom2, obs_domains.size());
  std::uint64_t tot_interactions{};
  std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom1, chrom2);
        const hictk::transformers::JoinGenomicCoords jsel(
            sel.template begin<std::uint64_t>(), sel.template end<std::uint64_t>(), fp.bins_ptr());

        for (const auto &p : jsel) {
          const auto delta = p.coords.bin2.start() - p.coords.bin1.start();

          const auto i1 = p.coords.bin1.rel_id();
          const auto i2 = p.coords.bin2.rel_id();

          if (delta < min_delta || delta >= max_delta || bin1_mask[i1] || bin2_mask[i2])
              [[unlikely]] {
            continue;
          }

          if (obs_domains.add_interactions(p) != 0) {
            tot_interactions += p.count;
            const auto exp_count = expected_matrix.at(i1, i2);
            exp_domains.add_interactions(hictk::Pixel{p.coords, exp_count});
          }
        }
      },
      f.get());

  return tot_interactions;
}

[[nodiscard]] static std::uint64_t map_interactions_to_domains_multi_pass(
    const hictk::File &f, GenomicDomainsIndexed<std::uint64_t> &obs_domains,
    GenomicDomainsIndexed<double> &exp_domains, const ExpectedMatrixStats &expected_matrix,
    std::uint64_t min_delta, std::uint64_t max_delta, const std::vector<bool> &bin1_mask,
    const std::vector<bool> &bin2_mask) {
  const auto &chrom1 = obs_domains.chrom1().name();
  const auto &chrom2 = obs_domains.chrom2().name();

  SPDLOG_DEBUG(
      "[{}:{}]: mapping interactions to {} genomic domains using the multi-pass strategy...",
      chrom1, chrom2, obs_domains.size());

  std::uint64_t tot_interactions{};
  const auto domains = obs_domains.to_vector();

  std::visit(
      [&](const auto &fp) {
        for (const auto &dom : domains | std::views::keys) {
          const auto sel =
              fp.fetch(chrom1, dom.start1(), dom.end1(), chrom2, dom.start2(), dom.end2());
          const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<std::uint64_t>(),
                                                            sel.template end<std::uint64_t>(),
                                                            fp.bins_ptr());

          for (const auto &p : jsel) {
            const auto delta = p.coords.bin2.start() - p.coords.bin1.start();

            const auto i1 = p.coords.bin1.rel_id();
            const auto i2 = p.coords.bin2.rel_id();

            if (delta < min_delta || delta >= max_delta || bin1_mask[i1] || bin2_mask[i2])
                [[unlikely]] {
              continue;
            }

            if (obs_domains.add_interactions(p) != 0) {
              tot_interactions += p.count;
              const auto exp_count = expected_matrix.at(i1, i2);
              exp_domains.add_interactions(hictk::Pixel{p.coords, exp_count});
            }
          }
        }
      },
      f.get());

  return tot_interactions;
}

[[nodiscard]] static std::vector<std::tuple<BEDPE, std::uint64_t, double>>
map_interactions_to_domains(const hictk::File &f, const GenomicDomains &domains,
                            const ExpectedMatrixStats &expected_matrix,
                            const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
                            std::uint64_t min_delta, std::uint64_t max_delta,
                            const std::vector<bool> &bin1_mask, const std::vector<bool> &bin2_mask,
                            DomainAggregationStrategy aggregation_stategy) {
  const auto t0 = std::chrono::system_clock::now();

  if (chrom1 != chrom2) {
    min_delta = 0;
    max_delta = std::numeric_limits<std::uint64_t>::max();
  }

  auto obs_domains = domains.fetch<std::uint64_t>(chrom1, chrom2);
  if (obs_domains.empty()) {
    return {};
  }

  auto exp_domains = domains.fetch<double>(chrom1, chrom2);

  if (aggregation_stategy == DomainAggregationStrategy::AUTO) {
    const auto coverage = obs_domains.coverage(std::max(std::uint32_t{100'000}, f.resolution()));
    SPDLOG_DEBUG("[{}:{}]: {} domain(s) cover ~{:.2f}% of the matrix", chrom1.name(), chrom2.name(),
                 obs_domains.size(), 100 * coverage);
    aggregation_stategy = coverage > 0.33 ? DomainAggregationStrategy::SINGLE_PASS
                                          : DomainAggregationStrategy::MULTI_PASS;
  }

  const auto tot_interactions =
      aggregation_stategy == DomainAggregationStrategy::SINGLE_PASS
          ? map_interactions_to_domains_one_pass(f, obs_domains, exp_domains, expected_matrix,
                                                 min_delta, max_delta, bin1_mask, bin2_mask)
          : map_interactions_to_domains_multi_pass(f, obs_domains, exp_domains, expected_matrix,
                                                   min_delta, max_delta, bin1_mask, bin2_mask);

  std::vector<std::tuple<BEDPE, std::uint64_t, double>> results(obs_domains.size());
  auto obs_domains_sorted = obs_domains.to_vector();
  const auto exp_domains_sorted = exp_domains.to_vector();
  assert(obs_domains_sorted.size() == exp_domains_sorted.size());

  for (std::size_t i = 0; i < obs_domains_sorted.size(); ++i) {
    results[i] = std::make_tuple(std::move(obs_domains_sorted[i].first),
                                 obs_domains_sorted[i].second, exp_domains_sorted[i].second);
  }

  [[maybe_unused]] const auto t1 = std::chrono::system_clock::now();
  SPDLOG_DEBUG("[{}:{}]: mapped {} interactions to {} genomic domains in {}!", chrom1.name(),
               chrom2.name(), tot_interactions, obs_domains.size(), format_duration(t1 - t0));

  return results;
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

[[nodiscard]] static std::size_t run_nchg_compute_worker(
    const ComputePvalConfig &c, const std::optional<GenomicDomains> &domains,
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

    const auto domains_ =
        parse_domains(f->chromosomes(), c.path_to_domains, true, true, chrom1, chrom2);
    return process_domains(f, domains_, expected_values, c);
  }

  return process_one_chromosome_pair(f, expected_values, c);
}

[[nodiscard]] static ChromosomePairs init_cis_chromosomes(
    const hictk::Reference &chroms, const std::optional<GenomicDomains> &domains) {
  ChromosomePairs buffer{};

  if (domains.has_value()) {
    for (const auto &domain : (*domains)()) {
      const auto &domain1 = domain.range1();
      const auto &domain2 = domain.range2();
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
    const hictk::Reference &chroms, const std::optional<GenomicDomains> &domains) {
  ChromosomePairs buffer{};

  if (domains.has_value()) {
    for (const auto &domain : (*domains)()) {
      const auto &domain1 = domain.range1();
      const auto &domain2 = domain.range2();
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

class MessageQueue {
  std::string _name;
  std::unique_ptr<boost::interprocess::message_queue> _queue{};
  std::mutex _mtx;
  std::string _eoq_signal{};
  const std::atomic<bool> *_early_return{};
  bool _destroy_queue{false};

  static constexpr std::size_t _max_message_size = sizeof(char) * 2048UL;

  MessageQueue(std::string name, std::size_t num_proc, const std::atomic<bool> *early_return,
               bool create)
      : _name(std::move(name)),
        _queue(create ? create_queue(_name, num_proc) : open_queue(_name)),
        _eoq_signal(
            fmt::format("### EOQ {}", std::chrono::system_clock::now().time_since_epoch().count())),
        _early_return(early_return),
        _destroy_queue(create) {}

 public:
  MessageQueue() = default;

  [[nodiscard]] static MessageQueue create(const std::string &name, std::size_t num_proc,
                                           const std::atomic<bool> &early_return) {
    return {name, num_proc, &early_return, true};
  }

  [[nodiscard]] static MessageQueue open(const std::string &name) {
    return {name, 1, nullptr, false};
  }

  MessageQueue(const MessageQueue &other) = delete;
  MessageQueue(MessageQueue &&other) noexcept
      : _name(std::move(other._name)),
        _queue(std::move(other._queue)),
        _eoq_signal(std::move(other._eoq_signal)),
        _early_return(other._early_return) {}

  ~MessageQueue() noexcept { close(); }

  MessageQueue &operator=(const MessageQueue &other) = delete;
  MessageQueue &operator=(MessageQueue &&other) noexcept {
    if (this == &other) {
      return *this;
    }

    _name = std::move(other._name);
    _queue = std::move(other._queue);
    _eoq_signal = std::move(other._eoq_signal);
    _early_return = other._early_return;

    return *this;
  }

  [[nodiscard]] std::string_view name() const noexcept { return _name; }

  void send(const spdlog::details::log_msg &msg) {
    assert(_queue);
    std::string buff{};
    std::unique_lock lck(_mtx);
    if (!glz::write_beve(msg, buff)) [[likely]] {
      send(buff, std::move(lck));
    }
  }

  void send(std::string_view msg) { send(msg, std::unique_lock{_mtx}); }

  void send(std::string_view msg, std::unique_lock<std::mutex> lck) {
    assert(_queue);
    assert(lck.owns_lock());

    while (true) {
      if (_early_return && *_early_return) {
        return;
      }

      if (!lck.owns_lock()) {
        lck.lock();
      }

      static const std::chrono::microseconds wait_time{50'000};
      static const boost::interprocess::ustime wait_time_us{
          static_cast<std::uint64_t>(wait_time.count())};

      if (!msg.empty() && msg.size() <= _max_message_size) {
        const auto sent = _queue->timed_send(msg.data(), msg.size(), 0, wait_time_us);
        if (sent) {
          return;
        }
        lck.unlock();
        std::this_thread::sleep_for(wait_time);
      }
    }
  }

  bool receive() {
    assert(_queue);
    assert(_early_return);

    std::string buff(_max_message_size, '\0');
    std::size_t msg_size{};
    [[maybe_unused]] std::uint32_t _{};

    std::unique_lock lck(_mtx, std::defer_lock);
    static const std::chrono::microseconds wait_time{50'000};
    static const boost::interprocess::ustime wait_time_us{
        static_cast<std::uint64_t>(wait_time.count())};

    while (true) {
      if (*_early_return) {
        return false;
      }

      lck.lock();
      const auto received =
          _queue->timed_receive(buff.data(), buff.size(), msg_size, _, wait_time_us);
      if (received) {
        break;
      }

      lck.unlock();
      std::this_thread::sleep_for(wait_time);
    }

    buff.resize(msg_size);
    if (buff == _eoq_signal) {
      return false;
    }

    spdlog::details::log_msg msg;
    if (glz::read_beve(msg, buff)) [[unlikely]] {
      fmt::println(stderr, "{}", buff);
    } else {
      spdlog::default_logger_raw()->log(
          msg.time, msg.source, msg.level != spdlog::level::info ? msg.level : spdlog::level::debug,
          msg.payload);
    }

    return true;
  }

  void send_eoq_signal() { send(_eoq_signal); }

  void close() noexcept {
    try {
      [[maybe_unused]] const std::scoped_lock lock(_mtx);
      _queue.reset();
      if (_destroy_queue) {
        assert(!_name.empty());
        boost::interprocess::message_queue::remove(_name.c_str());
      }
    } catch (...) {  // NOLINT
    }
  }

 private:
  [[nodiscard]] static std::unique_ptr<boost::interprocess::message_queue> open_queue(
      const std::filesystem::path &path) {
    SPDLOG_DEBUG("opening message queue with name={}...", path);
    return std::make_unique<boost::interprocess::message_queue>(boost::interprocess::open_only,
                                                                path.c_str());
  }
  [[nodiscard]] static std::unique_ptr<boost::interprocess::message_queue> create_queue(
      const std::filesystem::path &path, std::size_t num_proc) {
    assert(num_proc != 0);
    SPDLOG_DEBUG("initializing a message queue with name={}...", path);
    boost::interprocess::message_queue::remove(path.c_str());
    return std::make_unique<boost::interprocess::message_queue>(
        boost::interprocess::create_only, path.c_str(), num_proc * 4, sizeof(char) * 2048UL);
  }
};

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

static void process_log_messages(MessageQueue &msg_queue, std::atomic<bool> &early_return) {
  std::size_t num_except = 0;
  SPDLOG_DEBUG("starting logger thread...");

  for (std::size_t i = 0; !early_return; ++i) {
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
      const auto &[chrom1, chrom2] = chrom_pair;
      SPDLOG_DEBUG("submitting task {}/{} ({}:{})...", tasks_submitted + 1, workers.size(),
                   chrom1.name(), chrom2.name());
      return worker_fx(file_store, chrom1, chrom2, domains, tmpdir, msg_queue, ctx, base_config,
                       user_provided_expected_values, early_return);
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
                                   ? run_nchg_compute_worker(config, domains, expected_values)
                                   : run_nchg_compute_worker(config, domains);
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

[[nodiscard]] static ExpectedValues init_cis_expected_values(const ComputePvalConfig &c) {
  assert(c.compute_cis || c.chrom1 == c.chrom2);

  SPDLOG_INFO("initializing expected values for cis matrices...");
  const auto f = std::make_shared<hictk::File>(c.path_to_hic.string(), c.resolution);

  const auto bin_mask = parse_bin_mask(f->chromosomes(), f->resolution(), c.path_to_bin_mask);

  return ExpectedValues::cis_only(
      f,
      {.mad_max = c.mad_max,
       .min_delta = c.min_delta,
       .max_delta = c.max_delta,
       .bin_aggregation_possible_distances_cutoff = c.bin_aggregation_possible_distances_cutoff,
       .bin_aggregation_observed_distances_cutoff = c.bin_aggregation_observed_distances_cutoff,
       .interpolate = c.interpolate_expected_values,
       .interpolation_qtile = c.interpolation_qtile,
       .interpolation_window_size = c.interpolation_window_size},
      bin_mask);
}

static std::size_t process_queries(FileStore &file_store, const ChromosomePairs &chrom_pairs,
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

[[nodiscard]] static auto generate_execution_plan(const ComputePvalConfig &c,
                                                  bool init_file_store) {
  struct Plan {
    std::unique_ptr<FileStore> file_store;
    ChromosomePairs chrom_pairs{};
    std::optional<ExpectedValues> expected_values{};
    std::optional<GenomicDomains> domains{};
  };

  Plan plan{};

  const hictk::File f(c.path_to_hic.string(), c.resolution);

  if (f.bins().type() != hictk::BinTable::Type::fixed) {
    throw std::runtime_error("only file with uniform bin sizes are supported.");
  }

  if (!c.path_to_domains.empty()) {
    assert(!c.chrom1.has_value() || f.chromosomes().contains(*c.chrom1));
    assert(!c.chrom2.has_value() || f.chromosomes().contains(*c.chrom2));

    plan.domains.emplace(parse_domains(f.chromosomes(), c.path_to_domains, c.compute_cis,
                                       c.compute_trans, c.chrom1, c.chrom2));
  }

  if (!c.chrom1.has_value()) {
    assert(!c.chrom2.has_value());
    if (c.compute_cis) {
      plan.chrom_pairs = init_cis_chromosomes(f.chromosomes(), plan.domains);
    }

    if (c.compute_trans) {
      std::ranges::move(init_trans_chromosomes(f.chromosomes(), plan.domains),
                        std::inserter(plan.chrom_pairs, plan.chrom_pairs.end()));
    }
  } else {
    assert(c.chrom2.has_value());
    plan.chrom_pairs.emplace(
        std::make_pair(f.chromosomes().at(*c.chrom1), f.chromosomes().at(*c.chrom2)));
  }

  if (init_file_store) {
    process_file_collisions(c.output_prefix, plan.chrom_pairs, c.force);
  } else {
    assert(!c.output_path.empty());
    if (c.force) {
      const auto removed = std::filesystem::remove(c.output_path);  // NOLINT
      if (removed) {
        SPDLOG_DEBUG("file \"{}\" has been deleted", c.output_path.string());
      }
    } else if (std::filesystem::exists(c.output_path)) {
      throw std::runtime_error(fmt::format(
          "refusing to overwrite file \"{}\". Pass --force to overwrite.", c.output_path.string()));
    }
  }

  if (!c.path_to_expected_values.empty()) {
    assert(!plan.expected_values.has_value());
    plan.expected_values = ExpectedValues::deserialize(c.path_to_expected_values);
    validate_expected_values(*plan.expected_values, c.path_to_expected_values, plan.chrom_pairs,
                             f.resolution());
  }

  if (init_file_store) {
    const auto root_dir = c.output_prefix.has_parent_path() ? c.output_prefix.parent_path()
                                                            : std::filesystem::current_path();

    if (!root_dir.empty() && !std::filesystem::exists(root_dir)) {
      std::filesystem::create_directories(root_dir);  // NOLINT
    }
    plan.file_store = std::make_unique<FileStore>(
        root_dir, false, fmt::format("{}.json", c.output_prefix.filename().string()));
  }

  if (!plan.expected_values.has_value() && c.compute_cis) {
    plan.expected_values = init_cis_expected_values(c);
  }

  return plan;
}

template <typename... T>
static void warn_or_print(fmt::format_string<T...> fmt, T &&...args) {
  auto logger = spdlog::default_logger();
  if (logger) {
    SPDLOG_WARN(fmt, std::forward<T>(args)...);
  } else {
    fmt::print(stderr, fmt, std::forward<T>(args)...);
  }
}

static bool setup_file_backed_logger(const std::filesystem::path &output_prefix, bool force) {
  //                                  [2021-08-12 17:49:34.581] [info]: my log msg
  static constexpr auto *log_pattern{"[%Y-%m-%d %T.%e] %^[%l]%$: %v"};

  const auto log_file_path = fmt::format("{}.log", output_prefix.string());

  try {
    auto logger = spdlog::default_logger();
    if (!logger) {
      throw spdlog::spdlog_ex("application logger has not been properly setup");
    }

    if (!force && std::filesystem::exists(log_file_path)) {
      throw std::runtime_error("refusing to overwrite existing file: pass --force to overwrite.");
    }
    std::filesystem::remove(log_file_path);  // NOLINT

    if (output_prefix.has_parent_path()) {
      std::filesystem::create_directories(output_prefix.parent_path());  // NOLINT
    }

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file_path);
    file_sink->set_level(spdlog::level::debug);
    file_sink->set_pattern(log_pattern);

    auto tmp_logger = std::make_shared<spdlog::logger>("tmp_logger", file_sink);
    tmp_logger->set_level(spdlog::level::info);
    tmp_logger->log(spdlog::level::info, "log file was generated by NCHG v0.0.2");  // TODO fixme

    logger->sinks().emplace_back(std::move(file_sink));
    logger->set_level(spdlog::level::debug);

    return true;

  } catch (const spdlog::spdlog_ex &e) {
    warn_or_print("unable to initialize log file \"{}\": {}", log_file_path, e.what());
  } catch (const std::exception &e) {
    warn_or_print("unable to initialize log file \"{}\": {}", log_file_path, e.what());
  } catch (...) {
    warn_or_print("unable to initialize log file \"{}\": unknown error", log_file_path);
  }

  return false;
}

static void setup_beve_logger(const std::string &name) {
  auto logger = spdlog::default_logger();
  if (!logger || name.empty()) {
    return;
  }

  auto message_queue = std::make_shared<MessageQueue>(MessageQueue::open(name));

  auto callback_sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [message_queue](const spdlog::details::log_msg &msg) { message_queue->send(msg); });

  callback_sink->set_level(spdlog::level::debug);

  logger->sinks().clear();
  logger->sinks().emplace_back(std::move(callback_sink));
  logger->set_level(spdlog::level::debug);
}

int run_command(const ComputePvalConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();

  if (!c.log_message_queue.empty()) {
    setup_beve_logger(c.log_message_queue);
  }

  if (c.chrom1.has_value()) {
    assert(c.chrom2.has_value());
    assert(!c.output_path.empty());

    const auto plan = generate_execution_plan(c, false);

    const auto interactions_processed =
        run_nchg_compute_worker(c, plan.domains, plan.expected_values);
    const auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO("[{}:{}]: processed {} records in {}!", *c.chrom1, *c.chrom2,
                interactions_processed, format_duration(t1 - t0));

    SPDLOG_INFO("[{}:{}]: all records have been written to file \"{}\"", *c.chrom1, *c.chrom2,
                c.output_path.string());
    return 0;
  }

  const auto log_file_initialized = setup_file_backed_logger(c.output_prefix, c.force);

  auto [file_store, chrom_pairs, expected_values, domains] = generate_execution_plan(c, true);

  const auto interactions_processed =
      process_queries(*file_store, chrom_pairs, domains, expected_values, c);
  file_store->finalize();

  const auto t1 = std::chrono::steady_clock::now();
  if (interactions_processed == 0) {
    SPDLOG_WARN("no records have been processed. Is this intended?");
  } else {
    SPDLOG_INFO("processed {} records in {}!", interactions_processed, format_duration(t1 - t0));
  }

  auto files_created = chrom_pairs.size();
  ++files_created;  // chrom.sizes file
  ++files_created;  // report file
  files_created += static_cast<std::size_t>(log_file_initialized);

  SPDLOG_INFO("created {} new file(s) under prefix \"{}\"", files_created, c.output_prefix);

  return 0;
}

}  // namespace nchg
