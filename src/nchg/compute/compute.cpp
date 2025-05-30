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
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/bin_table.hpp>
#include <hictk/file.hpp>
#include <hictk/reference.hpp>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "./genomic_domains.hpp"
#include "./logging.hpp"
#include "./scheduler.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/file_store.hpp"
#include "nchg/text.hpp"
#include "nchg/tools/common.hpp"
#include "nchg/tools/config.hpp"

namespace nchg {

[[nodiscard]] static bool chrom_pair_is_empty(const hictk::File &f, const hictk::Chromosome &chrom1,
                                              const hictk::Chromosome &chrom2) {
  const auto sel = f.fetch(chrom1.name(), chrom2.name());
  return sel.begin<int>() == sel.end<int>();
}

template <typename Pairs>
[[nodiscard]] static ChromosomePairs drop_empty_pairs(const hictk::File &f, Pairs chrom_pairs) {
  // NOLINTNEXTLINE(*-const-correctness)
  std::vector chrom_pairs_flat(chrom_pairs.begin(), chrom_pairs.end());

  ChromosomePairs buff;
  bool cis = true;
  for (auto &&[chrom1, chrom2] : chrom_pairs) {
    cis &= chrom1 == chrom2;
    if (!chrom_pair_is_empty(f, chrom1, chrom2)) {
      buff.emplace(std::move(chrom1), std::move(chrom2));
    }
  }

  if (chrom_pairs_flat.size() != buff.size()) {
    SPDLOG_INFO(
        "dropped {}/{} {} chromosome pair(s) because their corresponding matrix has no "
        "interactions",
        chrom_pairs_flat.size() - buff.size(), chrom_pairs_flat.size(), cis ? "cis" : "trans");
  }

  return buff;
}

[[nodiscard]] static ChromosomePairs init_cis_chrom_pairs(const hictk::File &f,
                                                          const GenomicDomains &domains,
                                                          bool skip_empty_matrices) {
  ChromosomePairs buffer{};
  [[maybe_unused]] const auto &chroms = f.chromosomes();

  for (const auto &domain : domains()) {
    const auto &domain1 = domain.range1();
    const auto &domain2 = domain.range2();
    assert(chroms.contains(domain1.chrom().name()));
    assert(chroms.contains(domain2.chrom().name()));

    if (domain1.chrom() == domain2.chrom()) {
      buffer.emplace(domain1.chrom(), domain2.chrom());
    }
  }

  if (skip_empty_matrices) {
    return drop_empty_pairs(f, std::move(buffer));
  }

  return buffer;
}

[[nodiscard]] static ChromosomePairs init_cis_chrom_pairs(const hictk::File &f,
                                                          bool skip_empty_matrices) {
  ChromosomePairs buffer{};
  [[maybe_unused]] std::size_t num_pairs = 0;
  [[maybe_unused]] std::size_t dropped_pairs = 0;
  const auto &chroms = f.chromosomes();
  for (const auto &chrom : chroms) {
    if (chrom.is_all()) [[unlikely]] {
      continue;
    }

    const auto keep = !skip_empty_matrices || !chrom_pair_is_empty(f, chrom, chrom);
    num_pairs += keep;       // NOLINT(*-implicit-bool-conversion)
    dropped_pairs += !keep;  // NOLINT(*-implicit-bool-conversion)
    if (keep) {
      buffer.emplace(chrom, chrom);
    }
  }

  if (dropped_pairs != 0) {
    SPDLOG_INFO(
        "dropped {}/{} cis chromosome pair(s) because their corresponding matrix has no "
        "interactions",
        dropped_pairs, num_pairs + dropped_pairs);
  }

  return buffer;
}

[[nodiscard]] static ChromosomePairs init_trans_chrom_pairs(const hictk::File &f,
                                                            const GenomicDomains &domains,
                                                            bool skip_empty_matrices) {
  ChromosomePairs buffer{};
  [[maybe_unused]] const auto &chroms = f.chromosomes();
  for (const auto &domain : domains()) {
    const auto &domain1 = domain.range1();
    const auto &domain2 = domain.range2();
    assert(chroms.contains(domain1.chrom().name()));
    assert(chroms.contains(domain2.chrom().name()));

    if (domain1.chrom() != domain2.chrom()) {
      buffer.emplace(domain1.chrom(), domain2.chrom());
    }
  }

  if (skip_empty_matrices) {
    return drop_empty_pairs(f, std::move(buffer));
  }

  return buffer;
}

[[nodiscard]] static ChromosomePairs init_trans_chrom_pairs(const hictk::File &f,
                                                            bool skip_empty_matrices) {
  ChromosomePairs buffer{};
  const auto &chroms = f.chromosomes();
  [[maybe_unused]] std::size_t num_pairs = 0;
  [[maybe_unused]] std::size_t dropped_pairs = 0;
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
      const auto keep = !skip_empty_matrices || !chrom_pair_is_empty(f, chrom1, chrom2);
      num_pairs += keep;       // NOLINT(*-implicit-bool-conversion)
      dropped_pairs += !keep;  // NOLINT(*-implicit-bool-conversion)
      if (keep) {
        buffer.emplace(chrom1, chrom2);
      }
    }
  }

  if (dropped_pairs != 0) {
    SPDLOG_INFO(
        "dropped {}/{} trans chromosome pair(s) because their corresponding matrix has no "
        "interactions",
        dropped_pairs, num_pairs + dropped_pairs);
  }

  return buffer;
}

[[nodiscard]] static ChromosomePairs init_cis_chrom_pairs(
    const hictk::File &f, const std::optional<GenomicDomains> &domains, bool skip_empty_matrices) {
  auto pairs = domains.has_value() ? init_cis_chrom_pairs(f, *domains, skip_empty_matrices)
                                   : init_cis_chrom_pairs(f, skip_empty_matrices);
  if (pairs.empty()) {
    SPDLOG_WARN(
        "all cis chromosome pair(s) have been discarded because they have no interactions. "
        "Is this intended?");
  }
  return pairs;
}

[[nodiscard]] static ChromosomePairs init_trans_chrom_pairs(
    const hictk::File &f, const std::optional<GenomicDomains> &domains, bool skip_empty_matrices) {
  auto pairs = domains.has_value() ? init_trans_chrom_pairs(f, *domains, skip_empty_matrices)
                                   : init_trans_chrom_pairs(f, skip_empty_matrices);
  if (pairs.empty()) {
    SPDLOG_WARN(
        "all trans chromosome pair(s) have been discarded because they have no interactions. "
        "Is this intended?");
  }
  return pairs;
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

static void process_file_collision(const std::filesystem::path &path, bool force) {
  assert(!path.empty());
  if (force) {
    const auto removed = std::filesystem::remove(path);  // NOLINT
    if (removed) {
      SPDLOG_DEBUG("file \"{}\" has been deleted", path.string());
    }
  } else if (std::filesystem::exists(path)) {
    throw std::runtime_error(fmt::format(
        "refusing to overwrite file \"{}\". Pass --force to overwrite.", path.string()));
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

[[nodiscard]] static std::optional<GenomicDomains> try_parse_genomic_domains(
    const ComputePvalConfig &c, const hictk::Reference &chroms) {
  if (c.path_to_domains.empty()) {
    return {};
  }

  assert(!c.chrom1.has_value() || chroms.contains(*c.chrom1));
  assert(!c.chrom2.has_value() || chroms.contains(*c.chrom2));

  return parse_domains(chroms, c.path_to_domains, c.compute_cis, c.compute_trans, c.chrom1,
                       c.chrom2);
}

[[nodiscard]] static ChromosomePairs init_chromosome_pairs(
    const ComputePvalConfig &c, const hictk::File &f,
    const std::optional<GenomicDomains> &domains) {
  const auto &chroms = f.chromosomes();

  if (c.chrom1.has_value()) {
    assert(c.chrom2.has_value());
    const auto &chrom1 = chroms.at(*c.chrom1);
    const auto &chrom2 = chroms.at(*c.chrom2);

    if (chrom_pair_is_empty(f, chrom1, chrom2)) {
      return {};
    }
    return {std::make_pair(chrom1, chrom2)};
  }

  assert(!c.chrom2.has_value());
  ChromosomePairs chrom_pairs;

  if (c.compute_cis) {
    chrom_pairs = init_cis_chrom_pairs(f, domains, c.skip_empty_matrices);
  }

  if (c.compute_trans) {
    std::ranges::move(init_trans_chrom_pairs(f, domains, c.skip_empty_matrices),
                      std::inserter(chrom_pairs, chrom_pairs.end()));
  }

  return chrom_pairs;
}

[[nodiscard]] static std::optional<ExpectedValues> try_read_expected_values(
    const std::filesystem::path &path, const ChromosomePairs &chrom_pairs,
    std::uint32_t resolution) {
  if (path.empty()) {
    return {};
  }

  auto evs = ExpectedValues::deserialize(path);
  validate_expected_values(evs, path, chrom_pairs, resolution);
  return evs;
}

[[nodiscard]] static std::unique_ptr<FileStore> create_file_store(
    const std::filesystem::path &prefix) {
  const auto root_dir =
      prefix.has_parent_path() ? prefix.parent_path() : std::filesystem::current_path();

  if (!root_dir.empty() && !std::filesystem::exists(root_dir)) {
    std::filesystem::create_directories(root_dir);  // NOLINT
  }

  return std::make_unique<FileStore>(root_dir, false,
                                     fmt::format("{}.json", prefix.filename().string()));
}

[[nodiscard]] static bool plan_includes_cis_chrom_pairs(const ChromosomePairs &pairs) {
  for (const auto &[chrom1, chrom2] : pairs) {
    if (chrom1 == chrom2) {
      return true;
    }
  }

  return false;
}

[[nodiscard]] static auto generate_execution_plan(const ComputePvalConfig &c,
                                                  bool init_file_store) {
  struct Plan {
    std::unique_ptr<FileStore> file_store;
    ChromosomePairs chrom_pairs;
    std::optional<ExpectedValues> expected_values;
    std::optional<GenomicDomains> domains;
  };

  const hictk::File f(c.path_to_hic.string(), c.resolution);

  if (f.bins().type() != hictk::BinTable::Type::fixed) {
    throw std::runtime_error("only file with uniform bin sizes are supported.");
  }

  Plan plan{};
  plan.domains = try_parse_genomic_domains(c, f.chromosomes());
  plan.chrom_pairs = init_chromosome_pairs(c, f, plan.domains);

  if (init_file_store) {
    process_file_collisions(c.output_prefix, plan.chrom_pairs, c.force);
  } else {
    process_file_collision(c.output_path, c.force);
  }

  plan.expected_values =
      try_read_expected_values(c.path_to_expected_values, plan.chrom_pairs, f.resolution());

  if (init_file_store) {
    plan.file_store = create_file_store(c.output_prefix);
  }

  if (!plan.expected_values.has_value() && plan_includes_cis_chrom_pairs(plan.chrom_pairs)) {
    plan.expected_values = init_cis_expected_values(c);
  }

  return plan;
}

[[nodiscard]] static int process_one_chromosme_pair(
    const ComputePvalConfig &c, [[maybe_unused]] const std::chrono::steady_clock::time_point &t0) {
  assert(c.chrom1.has_value());
  assert(c.chrom2.has_value());
  assert(!c.output_path.empty());

  const auto plan = generate_execution_plan(c, false);

  if (c.skip_empty_matrices) {
    const hictk::File hf{c.path_to_hic.string(), c.resolution};
    const auto &chrom1 = hf.chromosomes().at(*c.chrom1);
    const auto &chrom2 = hf.chromosomes().at(*c.chrom2);
    if (chrom_pair_is_empty(hf, chrom1, chrom2)) {
      SPDLOG_WARN("[{}:{}]: matrix has no interactions: returning immediately!", *c.chrom1,
                  *c.chrom2);
      return 0;
    }
  }

  const auto interactions_processed =
      process_chromosome_pair(c, plan.domains, plan.expected_values);
  const auto t1 = std::chrono::steady_clock::now();
  SPDLOG_INFO("[{}:{}]: processed {} records in {}!", *c.chrom1, *c.chrom2, interactions_processed,
              format_duration(t1 - t0));

  SPDLOG_INFO("[{}:{}]: all records have been written to file \"{}\"", *c.chrom1, *c.chrom2,
              c.output_path.string());
  return 0;
}

[[nodiscard]] static int process_many_chromosome_pairs(
    const ComputePvalConfig &c, const std::chrono::steady_clock::time_point &t0) {
  const auto log_file_initialized = setup_file_backed_logger(c.output_prefix, c.force);

  auto [file_store, chrom_pairs, expected_values, domains] = generate_execution_plan(c, true);

  const auto interactions_processed =
      dispatch_queries(*file_store, chrom_pairs, domains, expected_values, c);
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

// NOLINTNEXTLINE(*-use-internal-linkage)
int run_command(const ComputePvalConfig &c) {
  const auto t0 = std::chrono::steady_clock::now();

  if (!c.log_message_queue.empty()) {
    setup_beve_logger(c.log_message_queue);
  }

  if (c.chrom1.has_value()) {
    return process_one_chromosme_pair(c, t0);
  }
  return process_many_chromosome_pairs(c, t0);
}

}  // namespace nchg
