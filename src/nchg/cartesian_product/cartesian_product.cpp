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

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
#include <iostream>
#include <istream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/hash.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

namespace nchg {

using Domain = decltype(hictk::GenomicInterval::parse_bed(std::declval<std::string>()));
using DomainPtr = std::shared_ptr<Domain>;

struct DomainHasher {
  using is_transparent = void;

  [[nodiscard]] static std::size_t operator()(const Domain& domain) noexcept {
    const auto& [chrom, start, end] = domain;
    return internal::hash_combine(0, chrom, start, end);
  }

  [[nodiscard]] static std::size_t operator()(const DomainPtr& ptr) noexcept {
    if (!ptr) [[unlikely]] {
      return std::hash<std::shared_ptr<Domain>>{}(ptr);
    }

    return operator()(*ptr);
  }
};

struct DomainEq {
  using is_transparent = void;

  [[nodiscard]] static bool operator()(const Domain& domain1, const Domain& domain2) noexcept {
    const auto& [chrom1, start1, end1] = domain1;
    const auto& [chrom2, start2, end2] = domain2;

    return chrom1 == chrom2 && start1 == start2 && end1 == end2;
  }

  [[nodiscard]] static bool operator()(const DomainPtr& ptr1, const Domain& domain2) noexcept {
    if (!ptr1) [[unlikely]] {
      return false;
    }

    return operator()(*ptr1, domain2);
  }

  [[nodiscard]] static bool operator()(const Domain& domain1, const DomainPtr& ptr2) noexcept {
    return operator()(ptr2, domain1);
  }

  [[nodiscard]] static bool operator()(const DomainPtr& ptr1, const DomainPtr& ptr2) noexcept {
    if (!ptr1 || !ptr2) [[unlikely]] {
      return ptr1 == ptr2;
    }

    return operator()(*ptr1, *ptr2);
  }
};

[[nodiscard]] static hictk::Reference parse_chrom_sizes(const std::filesystem::path& path,
                                                        const std::optional<std::string>& chrom1,
                                                        const std::optional<std::string>& chrom2) {
  if (chrom1.has_value()) {
    assert(chrom2.has_value());

    const std::vector names{*chrom1, *chrom2};
    const std::vector sizes(2, std::numeric_limits<std::uint32_t>::max());
    if (chrom1 == chrom2) {
      return {names.begin(), names.begin() + 1, sizes.begin()};
    }

    return {names.begin(), names.end(), sizes.begin()};
  }

  if (path.empty()) {
    return {};
  }

  auto chroms = hictk::Reference::from_chrom_sizes(path);
  if (!chroms.empty()) {
    return chroms;
  }

  throw std::runtime_error(fmt::format("unable to parse any chromosomes from file \"{}\"", path));
}

static void sort_domains(std::vector<Domain>& domains,
                         const std::vector<std::pair<std::string, std::size_t>>& chroms) {
  phmap::flat_hash_map<std::string, std::size_t> chrom_ranks(chroms.size());
  std::ranges::transform(
      chroms, std::inserter(chrom_ranks, chrom_ranks.end()),
      [&](const auto& kv) { return std::make_pair(kv.first, chrom_ranks.size()); });

  assert(!chroms.empty());
  std::ranges::sort(domains, [&](const Domain& lhs, const Domain& rhs) {
    const auto& i1 = chrom_ranks.at(std::get<0>(lhs));
    const auto& i2 = chrom_ranks.at(std::get<0>(rhs));
    if (i1 != i2) {
      return i1 < i2;
    }

    const auto start1 = std::get<1>(lhs);
    const auto start2 = std::get<1>(rhs);
    if (start1 != start2) {
      return start1 < start2;
    }

    const auto end1 = std::get<2>(lhs);
    const auto end2 = std::get<2>(rhs);
    return end1 < end2;
  });
}

static void sort_domains(std::vector<Domain>& domains, const hictk::Reference& chroms) {
  assert(!chroms.empty());
  std::ranges::sort(domains, [&](const Domain& lhs, const Domain& rhs) {
    const auto& chrom1 = chroms.at(std::get<0>(lhs));
    const auto& chrom2 = chroms.at(std::get<0>(rhs));
    if (chrom1 != chrom2) {
      return chrom1 < chrom2;
    }

    const auto start1 = std::get<1>(lhs);
    const auto start2 = std::get<1>(rhs);
    if (start1 != start2) {
      return start1 < start2;
    }

    const auto end1 = std::get<2>(lhs);
    const auto end2 = std::get<2>(rhs);
    return end1 < end2;
  });
}

static void detect_unordered_records(
    std::string_view chrom_name, std::size_t line_number,
    std::vector<std::pair<std::string, std::size_t>>& parsed_chromosomes) {
  if (parsed_chromosomes.empty()) [[unlikely]] {
    assert(line_number == 0);
    parsed_chromosomes.emplace_back(std::string{chrom_name}, line_number);
    return;
  }

  if (parsed_chromosomes.back().first == chrom_name) [[likely]] {
    return;
  }

  const auto it = std::ranges::find_if(parsed_chromosomes,
                                       [&](const auto& kv) { return kv.first == chrom_name; });

  if (it == parsed_chromosomes.end()) [[likely]] {
    parsed_chromosomes.emplace_back(std::string{chrom_name}, line_number);
    return;
  }

  std::vector<std::string> chrom_names{};
  std::transform(it, parsed_chromosomes.end(), std::back_inserter(chrom_names), [](const auto& kv) {
    return fmt::format("{} (first seen on line #{})", kv.first, kv.second);
  });
  throw std::runtime_error(
      fmt::format("domains are not sorted by their genomic coordinates: line #{} "
                  "references chromosome \"{}\", which was first encountered on line "
                  "\"{}\" and was followed by the chromosome(s) listed below:\n - {}",
                  line_number, it->first, it->second, fmt::join(chrom_names, "\n - ")));
}

enum class ParseStatus : std::uint_fast8_t { PARSED, SKIPPED, DUPLICATE };

[[nodiscard]] static ParseStatus parse_record(
    std::string_view line, std::size_t line_number, const hictk::Reference& reference,
    std::vector<std::pair<std::string, std::size_t>>& parsed_chromosomes,
    phmap::flat_hash_map<DomainPtr, std::size_t, DomainHasher, DomainEq>& domains) {
  auto record = hictk::GenomicInterval::parse_bed(line);
  if (!reference.empty() && !reference.contains(std::get<0>(record))) {
    SPDLOG_DEBUG("skipping line #{}", line_number);
    return ParseStatus::SKIPPED;
  }

  if (reference.empty()) {
    const auto& chrom = std::get<0>(record);
    detect_unordered_records(chrom, line_number, parsed_chromosomes);
  }

  if (domains.contains(record)) [[unlikely]] {
    SPDLOG_DEBUG("found duplicate record in line #{}", line_number);
    return ParseStatus::DUPLICATE;
  }

  domains.emplace(std::make_shared<Domain>(std::move(record)), domains.size());
  return ParseStatus::PARSED;
}

[[nodiscard]] static std::vector<Domain> parse_domains(std::istream& stream,
                                                       const hictk::Reference& reference) {
  // pair<chrom_name, line_of_first_encounter>
  std::vector<std::pair<std::string, std::size_t>> parsed_chromosomes{};
  phmap::flat_hash_map<DomainPtr, std::size_t, DomainHasher, DomainEq> unique_domains;

  std::string line;
  std::size_t duplicate_records{};
  std::size_t skipped_records{};
  std::size_t i = 0;

  try {
    assert(stream.good());
    for (; std::getline(stream, line); ++i) {
      if (line.empty()) [[unlikely]] {
        continue;
      }
      if (line.back() == '\r') [[unlikely]] {
        line.resize(line.size() - 1);
      }
      const auto status = parse_record(truncate_record<3>(line, '\t'), i, reference,
                                       parsed_chromosomes, unique_domains);
      switch (status) {
        using enum ParseStatus;
        case SKIPPED: {
          ++skipped_records;
          break;
        }
        case DUPLICATE: {
          ++duplicate_records;
          break;
        }
        case PARSED: {
          break;
        }
        default:
          unreachable_code();
      }
    }
  } catch (const std::exception& e) {
    if (!stream.eof()) {
      throw std::runtime_error(fmt::format("failed to parse line {}: {}", i, e.what()));
    }
  } catch (...) {
    throw std::runtime_error(fmt::format("failed to parse line {}: unknown error", i));
  }

  if (unique_domains.empty()) {
    throw std::runtime_error("unable to parse any domains");
  }

  if (duplicate_records != 0) {
    SPDLOG_WARN("parser discarded {} duplicate domains", duplicate_records);
  }

  if (skipped_records != 0) {
    assert(!reference.empty());
    SPDLOG_WARN(
        "parser skipped {} record(s) because they referred to chromosomes not listed in the given "
        "reference genome.",
        skipped_records);
  }

  std::vector<Domain> domains(unique_domains.size());

  for (auto& [domain, j] : unique_domains) {
    domains[j] = std::move(*domain);
  }

  SPDLOG_INFO("begin sorting domains...");
  if (!reference.empty()) {
    sort_domains(domains, reference);
  } else {
    sort_domains(domains, parsed_chromosomes);
  }
  SPDLOG_INFO("done sorting!");

  return domains;
}

struct ChromIndex {
  std::string chrom;
  std::size_t start_offset{};
  std::size_t end_offset{};
};

[[nodiscard]] static std::vector<ChromIndex> index_chromosomes(const std::vector<Domain>& domains) {
  std::vector<ChromIndex> index{};

  for (std::size_t i = 0; i < domains.size(); ++i) {
    const auto& chrom = std::get<0>(domains[i]);
    const auto i0 = i;
    auto i1 = i0 + 1;
    for (; i1 < domains.size(); ++i1) {
      if (chrom != std::get<0>(domains[i1])) {
        break;
      }
    }
#if defined(__apple_build_version__) && __apple_build_version__ < 16000000
    index.emplace_back(ChromIndex{chrom, i0, i1});
#else
    index.emplace_back(chrom, i0, i1);
#endif
    i = i1 - 1;
  }

  return index;
}

[[nodiscard]] static std::vector<Domain> parse_domains(const std::filesystem::path& path,
                                                       const hictk::Reference& chroms) {
  const auto read_from_stdin = path == "-" || path == "stdin";

  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);
  if (read_from_stdin) {
    SPDLOG_INFO("reading domains from stdin...");
  } else {
    SPDLOG_INFO("reading domains from file \"{}\"...", path.string());
    fs.open(path);
  }

  try {
    auto domains = fs.is_open() ? parse_domains(fs, chroms) : parse_domains(std::cin, chroms);

    if (read_from_stdin) {
      SPDLOG_INFO("successfully read {} domains from stdin!", domains.size());
    } else {
      SPDLOG_INFO("successfully read {} domains from file \"{}\"!", domains.size(), path.string());
    }

    return domains;
  } catch (const std::exception& e) {
    if (read_from_stdin) {
      throw std::runtime_error(fmt::format("failed to parse domains from stdin: {}", e.what()));
    }
    throw std::runtime_error(
        fmt::format("failed to parse domains from file \"{}\": {}", path.string(), e.what()));

  } catch (...) {
    if (read_from_stdin) {
      throw std::runtime_error(fmt::format("failed to parse domains from stdin: unknown error"));
    }
    throw std::runtime_error(
        fmt::format("failed to parse domains from file \"{}\": unknown error", path.string()));
  }
}

[[nodiscard]] static std::size_t print_chunk(const std::vector<Domain>& domains, std::size_t i0,
                                             std::size_t i1, std::size_t j0, std::size_t j1) {
  assert(i0 <= i1);
  assert(j0 <= j1);
  assert(j0 >= i0);

  std::size_t num_records = 0;
  for (std::size_t i = i0; i < i1; ++i) {
    const auto& [chrom1, start1, end1] = domains[i];
    for (std::size_t j = std::max(j0, i); j < j1; ++j) {
      const auto& [chrom2, start2, end2] = domains[j];
      fmt::print(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}\n"), chrom1, start1, end1, chrom2, start2,
                 end2);
      ++num_records;
    }
  }

  return num_records;
}

int run_command(const CartesianProductConfig& c) {
  assert(c.process_cis || c.process_trans);

  const auto t0 = std::chrono::steady_clock::now();

  const auto domains = parse_domains(c.path_to_domains,
                                     parse_chrom_sizes(c.path_to_chrom_sizes, c.chrom1, c.chrom2));

  const auto index = index_chromosomes(domains);

  std::size_t records_processed = 0;
  for (std::size_t i = 0; i < index.size(); ++i) {
    const auto& [chrom1, i0, i1] = index[i];
    if (c.chrom1.has_value() && *c.chrom1 != chrom1) {
      continue;
    }
    for (std::size_t j = i; j < index.size(); ++j) {
      const auto& [chrom2, j0, j1] = index[j];
      if (c.chrom2.has_value() && *c.chrom2 != chrom2) {
        continue;
      }

      if (!c.process_cis && chrom1 == chrom2) {
        SPDLOG_DEBUG("skipping domains for {}:{}...", chrom1, chrom2);
        continue;
      }
      if (!c.process_trans && chrom1 != chrom2) {
        SPDLOG_DEBUG("skipping domains for {}:{}...", chrom1, chrom2);
        continue;
      }

      SPDLOG_DEBUG("printing domains for {}:{}...", chrom1, chrom2);
      records_processed += print_chunk(domains, i0, i1, j0, j1);
    }
  }

  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("processed {} records in {}!", records_processed, format_duration(t1 - t0));

  return 0;
}

}  // namespace nchg
