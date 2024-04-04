// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#include "nchg/cli.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "nchg/config.hpp"

namespace nchg {

Cli::Cli(int argc, char **argv) : _argc(argc), _argv(argv), _exec_name(*argv) { make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return _subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(get_subcommand());
}

auto Cli::parse_arguments() -> Config {
  try {
    _cli.name(_exec_name);
    _cli.parse(_argc, _argv);

    if (_cli.get_subcommand("compute")->parsed()) {
      _subcommand = subcommand::compute;
    } else if (_cli.get_subcommand("filter")->parsed()) {
      _subcommand = subcommand::filter;
    } else {
      _subcommand = subcommand::help;
    }
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    _exit_code = _cli.exit(e);
    return _config;
  } catch (const std::exception &e) {
    _exit_code = 1;
    throw std::runtime_error(
        fmt::format(FMT_STRING("An unexpected error has occurred while parsing "
                               "CLI arguments: {}. If you see this "
                               "message, please file an issue on GitHub"),
                    e.what()));

  } catch (...) {
    _exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI "
        "arguments! If you see this message, please "
        "file an issue on GitHub");
  }
  validate_args();
  transform_args();

  _exit_code = 0;
  return _config;
}

int Cli::exit(const CLI::ParseError &e) const { return _cli.exit(e); }

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  switch (s) {
    case compute:
      return "compute";
    case filter:
      return "filter";
    default:
      assert(s == help);
      return "--help";
  }
}

void Cli::make_cli() {
  _cli.name(_exec_name);
  _cli.description("NCHG.");
  _cli.set_version_flag("-V,--version", "0.0.1");
  _cli.require_subcommand(1);

  make_compute_subcommand();
  make_filter_subcommand();
}

void Cli::make_compute_subcommand() {
  auto &sc =
      *_cli.add_subcommand(
               "compute", "Compute p-values for interactions from a .hic or .cool file using NCHG.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = ComputePvalConfig{};
           });

  _config = ComputePvalConfig{};
  auto &c = std::get<ComputePvalConfig>(_config);

  // clang-format off
  sc.add_option(
    "hic-matrix",
    c.path,
    "Path to a matrix in .hic, .cool or .mcool file with interactions to be processed.")
    ->check(IsValidHiCFile | IsValidCoolerFile | IsValidMultiresCoolerFile)
    ->required();
  sc.add_option(
    "--chrom1",
    c.chrom1,
    "Name of the first chromosome.\n"
    "Used to compute p-values only for a chromosome-chromosome matrix of interest.")
    ->capture_default_str();
  sc.add_option(
    "--chrom2",
    c.chrom2,
    "Name of the second chromosome.\n"
    "Used to compute p-values only for a chromosome-chromosome matrix of interest.")
    ->capture_default_str();
  sc.add_flag(
    "--cis-only",
    c.cis_only,
    "Process interactions from intra-chromosomal matrices only.")
    ->capture_default_str();
  sc.add_flag(
    "--trans-only",
    c.trans_only,
    "Process interactions from inter-chromosomal matrices only.")
    ->capture_default_str();
  sc.add_option(
    "--min-delta",
    c.min_delta,
    "Minimum distance from the diagonal required for interactions to be considered.")
    ->check(CLI::PositiveNumber)
    ->capture_default_str();
  sc.add_option(
    "--max-delta",
    c.max_delta,
    "Maximum distance from the diagonal required for interactions to be considered.")
    ->check(CLI::PositiveNumber)
    ->capture_default_str();
  sc.add_flag(
    "--write-header,!--no-write-header",
    c.write_header,
    "Write the file header to stdout.")
    ->capture_default_str();
  // clang-format on

  sc.get_option("--chrom2")->needs("--chrom1");

  sc.get_option("--cis-only")->excludes("--chrom1");
  sc.get_option("--cis-only")->excludes("--chrom2");
  sc.get_option("--cis-only")->excludes("--trans-only");

  sc.get_option("--trans-only")->excludes("--chrom1");
  sc.get_option("--trans-only")->excludes("--chrom2");

  _config = std::monostate{};
}

void Cli::make_filter_subcommand() {
  [[maybe_unused]] auto &sc =
      *_cli.add_subcommand("filter",
                           "Filter statistically significant interactions and correct "
                           "p-values for multiple hypothesis testing.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = FilterConfig{};
           });

  _config = FilterConfig{};
  [[maybe_unused]] auto &c = std::get<FilterConfig>(_config);

  // clang-format off
  sc.add_option(
      "tsv",
      c.path,
      "Path to the TSV produced by NCHG compute.")
      ->check(CLI::ExistingFile)
      ->required();
  sc.add_option(
      "--fdr",
       c.fdr,
      "FDR threshold used to identify significant interactions.")
      ->check(CLI::Bound(0.0, 1.0))
      ->capture_default_str();
  sc.add_option(
      "--log-ratio",
      c.log_ratio,
      "Log-ratio cutoff used to identify significant interactions.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  sc.add_flag(
      "--keep-non-significant,!--drop-non-significant",
      c.keep_non_significant,
      "Output non-significant interactions (i.e. ignore --fdr and --log-ratio cutoffs).")
      ->capture_default_str();
  sc.add_flag(
      "--write-header,!--no-write-header",
      c.write_header,
      "Write the file header to stdout.")
      ->capture_default_str();
  // clang-format on

  sc.get_option("--keep-non-significant")->excludes(sc.get_option("--fdr"));
  sc.get_option("--keep-non-significant")->excludes(sc.get_option("--log-ratio"));

  _config = std::monostate{};
}

void Cli::validate_args() const {
  switch (get_subcommand()) {
    case compute:
      return validate_compute_subcommand();  // NOLINT
    case filter:
      return validate_filter_subcommand();  // NOLINT
    case help:
      return;
  }
}

void Cli::validate_compute_subcommand() const {
  const auto &c = std::get<ComputePvalConfig>(_config);
  const auto &sc = *_cli.get_subcommand("compute");

  std::vector<std::string> warnings;
  std::vector<std::string> errors;

  if (c.min_delta >= c.max_delta) {
    errors.emplace_back("--min-delta should be less than --max-delta.");
  }

  const auto min_delta_parsed = !sc.get_option("--min-delta")->empty();
  const auto max_delta_parsed = !sc.get_option("--max-delta")->empty();
  const auto trans_only_parsed = !sc.get_option("--trans-only")->empty();

  if (trans_only_parsed && (min_delta_parsed || max_delta_parsed)) {
    warnings.emplace_back("--min-delta and --max-delta are ignored when --trans-only=true");
  }

  for (const auto &w : warnings) {
    SPDLOG_WARN(FMT_STRING("{}"), w);
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_filter_subcommand() const {}

void Cli::transform_args() {
  switch (get_subcommand()) {
    case compute:
      return transform_args_compute_subcommand();  // NOLINT
    case filter:
      return transform_args_filter_subcommand();  // NOLINT
    case help:
      return;
  }
}

void Cli::transform_args_compute_subcommand() {
  auto &c = std::get<ComputePvalConfig>(_config);
  if (c.chrom1 != "all" && c.chrom2 == "all") {
    c.chrom2 = c.chrom1;
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

void Cli::transform_args_filter_subcommand() {}

}  // namespace nchg
