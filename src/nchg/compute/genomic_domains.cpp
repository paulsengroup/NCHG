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

#include "./genomic_domains.hpp"

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/fmt.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/hash.hpp>
#include <hictk/pixel.hpp>
#include <hictk/reference.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/genomic_domains.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"

namespace nchg {

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

class IntervalStore {
  using GenomicInterval = std::shared_ptr<const hictk::GenomicInterval>;

  hictk::Reference _chroms;
  phmap::flat_hash_set<GenomicInterval, GIHasher, GIEqOperator> _intervals;

 public:
  IntervalStore() = default;
  explicit IntervalStore(hictk::Reference chroms) : _chroms(std::move(chroms)) {}

  [[nodiscard]] GenomicInterval emplace(std::shared_ptr<const hictk::GenomicInterval> domain) {
    assert(domain);
    assert(_chroms.contains(domain->chrom()));
    const auto match = _intervals.find(domain);
    if (match != _intervals.end()) {
      return *match;
    }

    const auto [it, _] = _intervals.emplace(std::move(domain));
    return *it;
  }

  [[nodiscard]] const hictk::Reference &chromosomes() const noexcept { return _chroms; }
};

class DomainParser {
  IntervalStore _intervals;
  std::optional<hictk::Chromosome> _chrom1;
  std::optional<hictk::Chromosome> _chrom2;

  bool _keep_cis{true};
  bool _keep_trans{true};

  std::size_t _domains_parsed{};
  std::size_t _domains_dropped{};
  std::size_t _dropped_cis{};
  std::size_t _dropped_trans{};

 public:
  DomainParser() = default;
  DomainParser(hictk::Reference chroms, std::optional<hictk::Chromosome> chrom1,
               std::optional<hictk::Chromosome> chrom2, bool keep_cis, bool keep_trans)
      : _intervals(std::move(chroms)),
        _chrom1(std::move(chrom1)),
        _chrom2(std::move(chrom2)),
        _keep_cis(keep_cis),
        _keep_trans(keep_trans) {}

  [[nodiscard]] std::optional<BEDPE> parse_domain(std::string_view buffer) {
    const auto record = truncate_record<6>(buffer);
    const auto domain1 = truncate_record<3>(record);
    const auto domain2 = truncate_record<3>(record.substr(domain1.size() + 1));

    const auto chrom1_parsed = truncate_record<1>(domain1);
    const auto chrom2_parsed = truncate_record<1>(domain2);

    ++_domains_parsed;

    if (!_keep_cis && chrom1_parsed == chrom2_parsed) {
      ++_dropped_cis;
      return {};
    }

    if (!_keep_trans && chrom1_parsed != chrom2_parsed) {
      ++_dropped_trans;
      return {};
    }

    const auto &chroms = _intervals.chromosomes();
    if (!_chrom1.has_value()) {
      assert(!_chrom2.has_value());
      if (!chroms.contains(chrom1_parsed) || !chroms.contains(chrom2_parsed)) {
        ++_domains_dropped;
        return {};
      }
    }

    if ((_chrom1.has_value() && *_chrom1 != chrom1_parsed) ||
        (_chrom2.has_value() && *_chrom2 != chrom2_parsed)) {
      ++_domains_dropped;
      return {};
    }

    auto gi1 = std::make_shared<const hictk::GenomicInterval>(
        hictk::GenomicInterval::parse_bed(chroms, domain1));
    auto gi2 = domain1 == domain2 ? gi1
                                  : std::make_shared<const hictk::GenomicInterval>(
                                        hictk::GenomicInterval::parse_bed(chroms, domain2));

    return BEDPE{_intervals.emplace(std::move(gi1)), _intervals.emplace(std::move(gi2))};
  }

  [[nodiscard]] constexpr std::size_t domains_parsed() const noexcept { return _domains_parsed; }
  [[nodiscard]] constexpr std::size_t dropped_cis() const noexcept { return _dropped_cis; }
  [[nodiscard]] constexpr std::size_t dropped_trans() const noexcept { return _dropped_trans; }
  [[nodiscard]] constexpr std::size_t domains_dropped() const noexcept { return _domains_dropped; }
  [[nodiscard]] std::optional<hictk::Chromosome> chrom1() const noexcept { return _chrom1; }
  [[nodiscard]] std::optional<hictk::Chromosome> chrom2() const noexcept { return _chrom2; }
};

static void generate_domain_parsing_report(const DomainParser &parser,
                                           std::size_t duplicate_domains) {
  if (duplicate_domains != 0) {
    SPDLOG_WARN("found {} duplicate domain(s)", duplicate_domains);
  }

  if (parser.domains_dropped() != 0) {
    if (parser.chrom1().has_value()) {
      // NOLINTBEGIN(*-unchecked-optional-access)
      assert(parser.chrom2().has_value());
      SPDLOG_DEBUG(
          "[{}:{}]: {}/{} domain(s) were dropped because they did not map to the specified "
          "chromosome(s)",
          parser.chrom1()->name(), parser.chrom2()->name(), parser.domains_dropped(),
          parser.domains_parsed());
      // NOLINTEND(*-unchecked-optional-access)
    } else {
      SPDLOG_WARN("{}/{} domain(s) were dropped because they did not map to any known chromosome",
                  parser.domains_dropped(), parser.domains_parsed());
    }
  }

  if (parser.dropped_cis() != 0) {
    SPDLOG_WARN(
        "{} domain(s) were dropped because they overlapped with the cis area of the interaction "
        "map",
        parser.dropped_cis());
  }
  if (parser.dropped_trans() != 0) {
    SPDLOG_WARN(
        "{} domain(s) were dropped because they overlapped with the trans area of the "
        "interaction map",
        parser.dropped_trans());
  }
}

GenomicDomains parse_domains(const hictk::Reference &chroms, const std::filesystem::path &path,
                             bool keep_cis, bool keep_trans,
                             const std::optional<hictk::Chromosome> &chrom1,
                             const std::optional<hictk::Chromosome> &chrom2) {
  if (chrom1.has_value()) {
    assert(chrom2.has_value());
    assert(chroms.contains(*chrom1));
    assert(chroms.contains(*chrom2));
  }

  SPDLOG_INFO("reading domains from file \"{}\"...", path.string());
  const auto t0 = std::chrono::steady_clock::now();

  std::ifstream ifs{};
  ifs.exceptions(ifs.exceptions() | std::ios::badbit | std::ios::failbit);

  std::size_t i = 1;
  std::size_t duplicate_domains{0};
  std::string buffer;

  DomainParser parser{chroms, chrom1, chrom2, keep_cis, keep_trans};
  phmap::flat_hash_set<BEDPE> domains_buffer{};

  try {
    ifs.open(path);

    for (; std::getline(ifs, buffer); ++i) {
      if (buffer.empty()) {
        continue;
      }

      if (buffer.back() == '\r') {
        buffer.resize(buffer.size() - 1);
      }

      auto domain = parser.parse_domain(buffer);
      if (!domain.has_value()) {
        continue;
      }

      if (domain->range1() > domain->range2()) {
        throw std::runtime_error(
            fmt::format("domains cannot overlap with the lower triangular matrix: offending "
                        "domain {:ucsc}; {:ucsc}",
                        domain->range1(), domain->range2()));
      }

      const auto &[_, inserted] = domains_buffer.emplace(std::move(*domain));
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

  generate_domain_parsing_report(parser, duplicate_domains);

  if (domains_buffer.empty()) {
    throw std::runtime_error(
        fmt::format("unable to parse any domain from file \"{}\"", path.string()));
  }

  GenomicDomains domains{std::vector<BEDPE>{std::make_move_iterator(domains_buffer.begin()),
                                            std::make_move_iterator(domains_buffer.end())},
                         true};

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("read {} domains from \"{}\" in {}", domains.size(), path.string(),
              format_duration(t1 - t0));
  return domains;
}

GenomicDomains parse_domains(const hictk::Reference &chroms, const std::filesystem::path &path,
                             bool keep_cis, bool keep_trans,
                             const std::optional<std::string> &chrom1,
                             const std::optional<std::string> &chrom2) {
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

std::filesystem::path write_domains_to_file(const GenomicDomains &domains,
                                            const std::filesystem::path &dest_dir,
                                            const hictk::Chromosome &chrom1,
                                            const hictk::Chromosome &chrom2, bool force) {
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

void write_chrom_sizes_to_file(const hictk::Reference &chroms, const std::filesystem::path &path,
                               bool force) {
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

[[nodiscard]] std::vector<std::tuple<BEDPE, std::uint64_t, double>> map_interactions_to_domains(
    const hictk::File &f, const GenomicDomains &domains, const ExpectedMatrixStats &expected_matrix,
    const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2, std::uint64_t min_delta,
    std::uint64_t max_delta, const std::vector<bool> &bin1_mask, const std::vector<bool> &bin2_mask,
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

}  // namespace nchg
