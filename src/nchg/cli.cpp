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

#include "nchg/tools/cli.hpp"

#include <arrow/io/file.h>
#include <fmt/format.h>
#include <parquet/arrow/reader.h>
#include <parquet/file_reader.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <exception>
#include <hictk/cooler/validation.hpp>
#include <hictk/hic/validation.hpp>
#include <hictk/numeric_utils.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>

#include <climits>

#else
#include <unistd.h>
#endif

#include "nchg/common.hpp"
#include "nchg/license.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/version.hpp"

namespace nchg {

template <typename N>
[[nodiscard]] static N parse_numeric(std::string_view tok) {
  return hictk::internal::parse_numeric_or_throw<N>(tok);
}

[[nodiscard]] static std::string pretty_format_memory(std::uint64_t mem, bool no_space = false) {
  // format memory with the appropriate unit and without uninformative decimal digits
  constexpr std::uint64_t kb = 1ULL << 10U;
  constexpr std::uint64_t mb = 1ULL << 20U;
  constexpr std::uint64_t gb = 1ULL << 30U;
  constexpr std::uint64_t tb = 1ULL << 40U;

  constexpr std::string_view kb_suffix_w_space = " KB";
  constexpr std::string_view mb_suffix_w_space = " MB";
  constexpr std::string_view gb_suffix_w_space = " GB";
  constexpr std::string_view tb_suffix_w_space = " TB";

  constexpr std::string_view kb_suffix_wo_space = "KB";
  constexpr std::string_view mb_suffix_wo_space = "MB";
  constexpr std::string_view gb_suffix_wo_space = "GB";
  constexpr std::string_view tb_suffix_wo_space = "TB";

  if (mem < kb) {
    return fmt::format("{} B", mem);
  }

  auto suffix = no_space ? tb_suffix_wo_space : tb_suffix_w_space;
  auto div = tb;

  if (mem < mb) {
    suffix = no_space ? kb_suffix_wo_space : kb_suffix_w_space;
    div = kb;
  } else if (mem < gb) {
    suffix = no_space ? mb_suffix_wo_space : mb_suffix_w_space;
    div = mb;
  } else if (mem < tb) {
    suffix = no_space ? gb_suffix_wo_space : gb_suffix_w_space;
    div = gb;
  }

  auto s = fmt::format("{:.3f}", static_cast<double>(mem) / static_cast<double>(div));
  while (!s.empty() && s.back() == '0') {
    s.resize(s.size() - 1);
  }

  if (!s.empty() && s.back() == '.') {
    s.resize(s.size() - 1);
  }

  if (s.empty()) [[unlikely]] {
    // this should never happen
    s = "0";
  }

  s.append(suffix);
  return s;
}

class MemoryValidator : public CLI::AsSizeValue {
  std::uint64_t _lb{};
  std::uint64_t _ub{std::numeric_limits<std::uint64_t>::max()};

 public:
  explicit MemoryValidator(std::uint64_t lb,
                           std::uint64_t ub = std::numeric_limits<std::uint64_t>::max())
      : AsSizeValue(false) {
    assert(lb <= ub);

    _lb = lb;
    _ub = ub;

    auto as_size_value_func = func_;

    func_ = [this, as_size_value_func](std::string &input) -> std::string {
      const auto original_input = input;
      if (const auto res = as_size_value_func(input); !res.empty()) {
        return res;
      }

      try {
        const auto mem = parse_numeric<std::uint64_t>(input);
        if (mem < _lb) {
          input = original_input;
          throw CLI::ValidationError(fmt::format("Memory cannot be less than {} (found {})",
                                                 pretty_format_memory(_lb),
                                                 pretty_format_memory(mem)));
        }

        if (mem > _ub) {
          input = original_input;
          throw CLI::ValidationError(fmt::format("Memory cannot be more than {} (found {})",
                                                 pretty_format_memory(_ub),
                                                 pretty_format_memory(mem)));
        }

        input = fmt::to_string(mem);
        return {};
      } catch (const CLI::ValidationError &) {
        throw;
      } catch (const std::exception &e) {
        throw CLI::ValidationError(
            fmt::format("Value \"{}\" is not a valid number: {}", input, e.what()));
      }
    };
  }

  [[nodiscard]] std::uint64_t lb() const noexcept { return _lb; }
  [[nodiscard]] std::uint64_t ub() const noexcept { return _ub; }
};

class ParquetFileValidator : public CLI::detail::ExistingFileValidator {
 public:
  ParquetFileValidator() {
    auto existing_file_func = func_;
    func_ = [existing_file_func](std::string &input) -> std::string {
      try {
        if (const auto res = existing_file_func(input); !res.empty()) {
          return res;
        }
        std::shared_ptr<arrow::io::ReadableFile> fp;
        PARQUET_ASSIGN_OR_THROW(fp, arrow::io::ReadableFile::Open(input))
        const auto result = parquet::arrow::OpenFile(fp, arrow::default_memory_pool());
        if (!result.ok()) {
          throw std::runtime_error(result.status().ToString());
        }
        return {};
      } catch (const std::exception &e) {
        throw CLI::ValidationError(
            fmt::format("failed to open file \"{}\" for reading: {}", input, e.what()));
      }
    };
  }
};

// NOLINTBEGIN(cert-err58-cpp)
// https://duckdb.org/docs/stable/guides/performance/environment.html#minimum-required-memory
static const auto IsValidMemoryDuckDB = MemoryValidator(125ULL << 20U);
static const auto IsValidParquet = ParquetFileValidator();
// NOLINTEND(cert-err58-cpp)

Cli::Cli(int argc, char **argv) : _argc(argc), _argv(argv), _exec_name(*argv) { make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return _subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(get_subcommand());
}

auto Cli::parse_arguments() -> Config {
  try {
    _cli.name(_exec_name);
    _cli.parse(_argc, _argv);

    if (handle_help_flags()) {
      return _config;
    }

    using enum subcommand;
    if (_cli.get_subcommand("cartesian-product")->parsed()) {
      _subcommand = cartesian_product;
    } else if (_cli.get_subcommand("checksum")->parsed()) {
      _subcommand = checksum;
    } else if (_cli.get_subcommand("compute")->parsed()) {
      _subcommand = compute;
    } else if (_cli.get_subcommand("expected")->parsed()) {
      _subcommand = expected;
    } else if (_cli.get_subcommand("filter")->parsed()) {
      _subcommand = filter;
    } else if (_cli.get_subcommand("merge")->parsed()) {
      _subcommand = merge;
    } else if (_cli.get_subcommand("metadata")->parsed()) {
      _subcommand = metadata;
    } else if (_cli.get_subcommand("view")->parsed()) {
      _subcommand = view;
    } else {
      _subcommand = none;
      for (const auto *opt : {"--help", "--version"}) {
        if (!_cli.get_option(opt)->empty()) {
          _exit_code = 0;
          return _config;
        }
      }
      fmt::print(stderr, "A subcommand is required\nRun with --help for more information.\n");
      _exit_code = 1;
      return _config;
    }
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    _exit_code = _cli.exit(e);
    return _config;
  } catch (const std::exception &e) {
    _exit_code = 1;
    throw std::runtime_error(
        fmt::format("An unexpected error has occurred while parsing "
                    "CLI arguments: {}. If you see this "
                    "message, please file an issue on GitHub",
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
int Cli::exit() const noexcept { return _exit_code; }

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  using enum subcommand;
  switch (s) {
    case cartesian_product:
      return "cartesian-product";
    case checksum:
      return "checksum";
    case compute:
      return "compute";
    case expected:
      return "expected";
    case filter:
      return "filter";
    case merge:
      return "merge";
    case metadata:
      return "metadata";
    case view:
      return "view";
    default:
      assert(s == none);
      return "";
  }
}

void Cli::log_warnings() const noexcept {
  for (const auto &w : _warnings) {
    SPDLOG_WARN("{}", w);
  }
  _warnings.clear();
}

void Cli::make_cli() {
  _cli.name(_exec_name);
  _cli.description("NCHG.");
  _cli.set_version_flag("-V,--version", std::string{config::version::str_long()});

  auto *grp = _cli.add_option_group("help");
  grp->require_option(0, 1);
  grp->set_help_flag();

  /*
  grp->add_flag_callback(
      "--help-cite", [this]() { _help_flag = "cite"; },
      "Print NCHG's citation in Bibtex format and exit.");
  */
  grp->add_flag_callback(
      "--help-docs", [this]() { _help_flag = "docs"; },
      "Print the URL to NCHG's documentation and exit.");
  grp->add_flag_callback(
      "--help-license", [this]() { _help_flag = "license"; }, "Print the NCHG license and exit.");

  make_cartesian_product_subcommand();
  make_checksum_subcommand();
  make_compute_subcommand();
  make_expected_subcommand();
  make_filter_subcommand();
  make_merge_subcommand();
  make_metadata_subcommand();
  make_view_subcommand();
}

void Cli::make_cartesian_product_subcommand() {
  auto &sc =
      *_cli.add_subcommand("cartesian-product",
                           "Compute the cartesian product of domains from a list in BED3 format.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = CartesianProductConfig{};
           });

  _config = CartesianProductConfig{};
  auto &c = std::get<CartesianProductConfig>(_config);

  // clang-format off
  sc.add_option(
    "domains",
    c.path_to_domains,
    "Path to a BED3+ file with the list of domains to be processed.\n"
    "Domains should be sorted by chromosome when --chrom-sizes is not provided.\n"
    "Pass \"-\" or \"stdin\" if the domains should be read from stdin.")
    ->check(CLI::ExistingFile | CLI::IsMember{{"-", "stdin"}})
    ->required();
  sc.add_option(
    "-c,--chrom-sizes",
    c.path_to_chrom_sizes,
    "Path to .chrom.sizes file.\n"
    "Chromosomes will be used to sort domains prior to processing.")
    ->check(CLI::ExistingFile);
  sc.add_option(
    "--chrom1",
    c.chrom1,
    "Name of the first chromosome used to filter domains.")
    ->capture_default_str();
  sc.add_option(
    "--chrom2",
    c.chrom2,
    "Name of the second chromosome used to filter domains.")
    ->capture_default_str();
  sc.add_flag_function(
    "--cis-only",
    [&c](auto n) { if (n != 0) {c.process_trans = false;} },
    "Only output pair of domains corresponding to regions interacting in cis.")
    ->capture_default_str();
  sc.add_flag_function(
    "--trans-only",
    [&c](auto n) { if (n != 0) {c.process_cis = false; }},
    "Only output pair of domains corresponding to regions interacting in trans.")
    ->capture_default_str();
  sc.add_option(
    "-v,--verbosity",
    c.verbosity,
    "Set verbosity of output to the console.")
    ->check(CLI::Range(1, 5))
    ->capture_default_str();
  // clang-format on

  sc.get_option("--chrom-sizes")->excludes("--chrom1");
  sc.get_option("--chrom-sizes")->excludes("--chrom2");

  sc.get_option("--cis-only")->excludes("--chrom1");
  sc.get_option("--cis-only")->excludes("--chrom2");
  sc.get_option("--cis-only")->excludes("--trans-only");

  sc.get_option("--trans-only")->excludes("--chrom1");
  sc.get_option("--trans-only")->excludes("--chrom2");

  _config = std::monostate{};
}

void Cli::make_checksum_subcommand() {
  auto &sc =
      *_cli.add_subcommand("checksum", "Checksum one or more files produced by NCHG compute.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = ChecksumConfig{};
           });

  _config = ChecksumConfig{};
  auto &c = std::get<ChecksumConfig>(_config);

  // clang-format off
  sc.add_option(
      "files",
      c.files,
      "Path to one or more file to be checksummed.")
      ->check(CLI::ExistingFile)
      ->required();
  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 5))
      ->capture_default_str();
  //clang-format on

  _config = std::monostate{};
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
    "hic-file",
    c.path_to_hic,
    "Path to a matrix in .hic, .cool or .mcool file with the interactions to be processed.")
    ->check(IsValidHiCFile | IsValidCoolerFile | IsValidMultiresCoolerFile)
    ->required();
  sc.add_option(
    "output-prefix",
    c.output_prefix,
    "Path prefix to use for output.\n"
    "Depending on the parameters used to invoke NCHG, this will result in one or more\n"
    "files named like myprefix.chrA.chrB.parquet plus a file named myprefix.chrom.sizes.\n"
    "When --chrom1 and/or --chrom2 have been specified a single file named myprefix.parquet will be created.")
    ->required();
  sc.add_option(
    "--resolution",
    c.resolution,
    "Matrix resolution. Required when the input file is in .hic or .mcool format.");
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
  sc.add_option(
    "--expected-values",
    c.path_to_expected_values,
    "Path to a .h5 file with the expected values generated by NCHG expected.")
    ->check(CLI::ExistingFile);
  sc.add_option(
    "--domains",
    c.path_to_domains,
    "Path to a BEDPE file with the list of 2D domains to be processed.\n"
    "For each domain, NCHG will first fetch and aggregate interactions overlapping with the given coordinates.\n"
    "Then, NCHG will asses the statistical significance of the observed interactions after aggregation.")
    ->check(CLI::ExistingFile);
  sc.add_flag(
    "--skip-empty-matrices,!--keep-empty-matrices",
    c.skip_empty_matrices,
    "Control whether NCHG should create empty .parquet file(s) for chromosome-chromosome\n"
    "matrices with no interactions.")
    ->capture_default_str();
  sc.add_flag_function(
    "--cis-only",
    [&c](auto n) { if (n != 0) {c.compute_trans = false;} },
    "Process interactions from intra-chromosomal matrices only.")
    ->capture_default_str();
  sc.add_flag_function(
    "--trans-only",
    [&c](auto n) { if (n != 0) {c.compute_cis = false; }},
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
  sc.add_option(
    "--bin-aggregation-possible-distances-cutoff",
    c.bin_aggregation_possible_distances_cutoff,
    "Cutoff on the number of possible bin-pair combinations used to determine\n"
    "when to aggregate bins from the expected matrix.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_option(
    "--bin-aggregation-observed-distances-cutoff",
    c.bin_aggregation_observed_distances_cutoff,
    "Cutoff on the number of observed bin-pair combinations used to determine\n"
    "when to aggregate bins from the expected matrix.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_flag(
    "--interpolate-expected-values,!--no-interpolate-expected-values",
    c.interpolate_expected_values,
    "Interpolate expected values profile to deal with outliers due to e.g. small translocations.")
    ->capture_default_str();
  sc.add_option(
    "--evs-interpolation-qtile",
    c.interpolation_qtile,
    "Percentile used to detect outliers during interpolation.")
    ->check(CLI::Bound(0.0, 1.0))
    ->capture_default_str();
  sc.add_option(
    "--evs-interpolation-window",
    c.interpolation_window_size,
    "Window size in bps used to interpolate expected value profiles.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_option(
    "--mad-max",
    c.mad_max,
    "Cutoff used by the MAD-max filter to mask bad bins.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_option(
    "--bad-bin-fraction",
    c.bad_bin_fraction,
    "Largest fraction of masked bins required for a domain to be processed.\n"
    "Has no effect when --domains is not used.")
  ->check(CLI::Bound(0.0, 1.0))
  ->capture_default_str();
  sc.add_option(
    "--interaction-aggregation-strategy",
    c.domain_aggregation_stategy,
    "Strategy used to aggregate interactions using the provided genomic domains.\n"
    "Ignored when --domains has not been specified.\n"
    "The \"one-pass\" algorithm is optimized for scenarios where the provided domains cover a\n"
    "significant fraction of individual chromosome-chromosome matrix (e.g. 30% or more),\n"
    "while the \"multi-pass\" algorithm is better at handling all other scenarios.")
  ->transform(
    CLI::CheckedTransformer{
      std::map<std::string, ComputePvalConfig::DomainAggregationStrategy>{
        {"auto", ComputePvalConfig::DomainAggregationStrategy::AUTO},
        {"one-pass", ComputePvalConfig::DomainAggregationStrategy::SINGLE_PASS},
        {"multi-pass", ComputePvalConfig::DomainAggregationStrategy::MULTI_PASS},
      }
    })
  ->default_str("auto");
  sc.add_option(
    "--bin-mask",
    c.path_to_bin_mask,
    "Path to a BED3+ file with a list of domains that should be masked.")
    ->check(CLI::ExistingFile)
    ->capture_default_str();
  sc.add_option(
    "--threads",
    c.threads,
    "Number of worker threads.")
    ->check(CLI::PositiveNumber)
    ->capture_default_str();
  sc.add_option(
    "--compression-level",
    c.compression_lvl,
    "Compression level used to compress columns in the output .parquet file.")
    ->check(CLI::Bound(1, 22))
    ->capture_default_str();
  sc.add_option(
    "--compression-method",
    c.compression_method,
    "Method used to compress individual columns in the .parquet file.")
    ->check(CLI::IsMember({"zstd", "lz4"}))
    ->capture_default_str();
  sc.add_flag(
    "--force",
    c.force,
    "Force overwrite existing output file(s).")
    ->capture_default_str();
  sc.add_option(
    "-v,--verbosity",
    c.verbosity,
    "Set verbosity of output to the console.")
    ->check(CLI::Range(1, 5))
    ->capture_default_str();
  // clang-format on

  sc.get_option("--chrom2")->needs("--chrom1");

  sc.get_option("--cis-only")->excludes("--chrom1");
  sc.get_option("--cis-only")->excludes("--chrom2");
  sc.get_option("--cis-only")->excludes("--trans-only");

  sc.get_option("--trans-only")->excludes("--chrom1");
  sc.get_option("--trans-only")->excludes("--chrom2");

  sc.get_option("--expected-values")->excludes("--min-delta");
  sc.get_option("--expected-values")->excludes("--max-delta");
  sc.get_option("--expected-values")->excludes("--mad-max");
  sc.get_option("--expected-values")->excludes("--bin-aggregation-possible-distances-cutoff");
  sc.get_option("--expected-values")->excludes("--bin-aggregation-observed-distances-cutoff");
  sc.get_option("--expected-values")->excludes("--interpolate-expected-values");
  sc.get_option("--expected-values")->excludes("--evs-interpolation-qtile");
  sc.get_option("--expected-values")->excludes("--evs-interpolation-window");
  sc.get_option("--expected-values")->excludes("--bin-mask");

  _config = std::monostate{};
}

void Cli::make_expected_subcommand() {
  [[maybe_unused]] auto &sc =
      *_cli.add_subcommand("expected", "Compute the expected profile for a given Hi-C matrix.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = ExpectedConfig{};
           });

  _config = ExpectedConfig{};
  auto &c = std::get<ExpectedConfig>(_config);

  // clang-format off
  sc.add_option(
    "hic-matrix",
    c.input_path,
    "Path to a matrix in .hic, .cool or .mcool file with interactions to be processed.")
    ->check(IsValidHiCFile | IsValidCoolerFile | IsValidMultiresCoolerFile)
    ->required();
  sc.add_option(
    "--output",
    c.output_path,
    "Path where to store the expected profiles.")
    ->required();
  sc.add_option(
    "--resolution",
    c.resolution,
    "Matrix resolution. Required when the input file is in .hic or .mcool format.");
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
    "Compute expected values from cis interactions only.")
    ->capture_default_str();
  sc.add_flag(
    "--trans-only",
    c.trans_only,
    "Compute expected values from trans interactions only.")
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
  sc.add_option(
    "--bin-aggregation-possible-distances-cutoff",
    c.bin_aggregation_possible_distances_cutoff,
    "Cutoff on the number of possible bin-pair combinations used to determine\n"
    "when to aggregate bins from the expected matrix.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_option(
    "--bin-aggregation-observed-distances-cutoff",
    c.bin_aggregation_observed_distances_cutoff,
    "Cutoff on the number of observed bin-pair combinations used to determine\n"
    "when to aggregate bins from the expected matrix.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_flag(
    "--interpolate-expected-values,!--no-interpolate-expected-values",
    c.interpolate_expected_values,
    "Interpolate expected values profile to deal with outliers due to e.g. small translocations.")
    ->capture_default_str();
  sc.add_option(
    "--evs-interpolation-qtile",
    c.interpolation_qtile,
    "Percentile used to detect outliers during interpolation.")
    ->check(CLI::Bound(0.0, 1.0))
    ->capture_default_str();
  sc.add_option(
    "--evs-interpolation-window",
    c.interpolation_window_size,
    "Window size in bps used to interpolate expected value profiles.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_option(
    "--mad-max",
    c.mad_max,
    "Cutoff used by the MAD-max filter to mask bad bins.")
    ->check(CLI::NonNegativeNumber)
    ->capture_default_str();
  sc.add_option(
    "--bin-mask",
    c.path_to_bin_mask,
    "Path to a BED3+ file with a list of domains that should be masked.")
    ->check(CLI::ExistingFile)
    ->capture_default_str();
  sc.add_flag(
    "--force",
    c.force,
    "Force overwrite existing output file(s).")
    ->capture_default_str();
  sc.add_option(
    "-v,--verbosity",
    c.verbosity,
    "Set verbosity of output to the console.")
    ->check(CLI::Range(1, 5))
    ->capture_default_str();
  // clang-format on

  sc.get_option("--chrom2")->needs("--chrom1");

  sc.get_option("--cis-only")->excludes("--trans-only");
  sc.get_option("--cis-only")->excludes("--chrom1");
  sc.get_option("--cis-only")->excludes("--chrom2");

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
  auto &c = std::get<FilterConfig>(_config);

  // clang-format off
  sc.add_option(
    "input-parquet",
    c.input_path,
    "Path to a parquet file produced by NCHG merge or compute.")
    ->check(IsValidParquet)
    ->required();
  sc.add_option(
    "output-parquet",
    c.output_path,
    "Path where to store the output table in parquet format.")
    ->required();
  sc.add_flag(
    "--force",
    c.force,
    "Force overwrite existing output file(s).")
    ->capture_default_str();
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
    "--drop-non-significant,!--keep-non-significant",
    c.drop_non_significant,
    "Drop non-significant interactions based on the --fdr and --log-ratio cutoffs.")
    ->capture_default_str();
  sc.add_flag(
    "--correct-pvals-chrom-by-chrom",
    c.correct_chrom_chrom_separately,
    "Perform multiple hypothesis correction by treating each pair of chromosomes as independent experiments.")
    ->capture_default_str();
  sc.add_flag(
    "--correct-pvals-cis-trans,!--correct-pvals-all",
    c.correct_cis_trans_separately,
    "Perform multiple hypothesis correction by treating cis and trans matrices as two separate experiments.")
    ->capture_default_str();
  sc.add_option(
    "--compression-level",
    c.compression_lvl,
    "Compression level used to compress columns in the output .parquet file.")
    ->check(CLI::Bound(1, 22))
    ->capture_default_str();
  sc.add_option(
    "--compression-method",
    c.compression_method,
    "Method used to compress individual columns in the .parquet file.")
    ->check(CLI::IsMember({"zstd", "lz4"}))
    ->capture_default_str();
  sc.add_option(
    "--threads",
    c.threads,
    "Number of worker threads used for compression.")
    ->check(CLI::Bound(2U, std::thread::hardware_concurrency()))
    ->capture_default_str();
  sc.add_option(
    "-v,--verbosity",
    c.verbosity,
    "Set verbosity of output to the console.")
    ->check(CLI::Range(1, 5))
    ->capture_default_str();
  // clang-format on

  sc.get_option("--keep-non-significant")->excludes(sc.get_option("--fdr"));
  sc.get_option("--keep-non-significant")->excludes(sc.get_option("--log-ratio"));
  sc.get_option("--correct-pvals-chrom-by-chrom")
      ->excludes(sc.get_option("--correct-pvals-cis-trans"));

  _config = std::monostate{};
}

void Cli::make_merge_subcommand() {
  [[maybe_unused]] auto &sc =
      *_cli.add_subcommand("merge", "Merge the output of NCHG compute into a single parquet file.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = MergeConfig{};
           });

  _config = MergeConfig{};
  auto &c = std::get<MergeConfig>(_config);

  auto *inputs_grp = sc.add_option_group("inputs");
  inputs_grp->require_option(1);
  auto *params_grp = sc.add_option_group("parameters");
  inputs_grp->set_help_flag();
  params_grp->set_help_flag();

  // clang-format off
  inputs_grp->add_option(
    "--input-prefix",
    c.input_prefix,
    "Path prefix where the files produced by NCHG compute are located.");
  inputs_grp->add_option(
    "--input-files",
    c.input_files,
    "Two or more paths to .parquet files generated by NCHG.")
    ->check(IsValidParquet);
  sc.add_option(
    "-o,--output",
    c.output_path,
    "Output path.")
    ->required();
  params_grp->add_flag(
    "--force",
    c.force,
    "Force overwrite existing output file(s).")
    ->capture_default_str();
  params_grp->add_flag(
    "--ignore-report-file,!--use-report-file",
    c.ignore_report_file,
    "Control whether the report file generated by NCHG compute should be\n"
    "used to validate the input files before merging.")
    ->capture_default_str();
  params_grp->add_option(
    "--compression-level",
    c.compression_lvl,
    "Compression level used to compress columns in the output .parquet file.")
    ->check(CLI::Bound(1, 22))
    ->capture_default_str();
  params_grp->add_option(
    "--compression-method",
    c.compression_method,
    "Method used to compress individual columns in the .parquet file.")
    ->check(CLI::IsMember({"zstd", "lz4"}))
    ->capture_default_str();
  params_grp->add_option(
    "--threads",
    c.threads,
    "Number of worker threads.")
    ->check(CLI::Range(2U, std::max(2U, std::thread::hardware_concurrency())))
    ->capture_default_str();
  params_grp->add_option(
    "--max-memory-per-thread",
    c.memory_per_thread,
    "Memory budget for each processing thread.")
    ->transform(IsValidMemoryDuckDB)  // NOLINT(cppcoreguidelines-slicing)
    ->default_str(pretty_format_memory(c.memory_per_thread, true));
  params_grp->add_option(
    "--max-memory",
    c.memory_limit,
    "Memory budget for the whole process.")
    ->transform(IsValidMemoryDuckDB)  // NOLINT(cppcoreguidelines-slicing)
    ->capture_default_str();
  params_grp->add_option(
    "-v,--verbosity",
    c.verbosity,
    "Set verbosity of output to the console.")
    ->check(CLI::Range(1, 5))
    ->capture_default_str();
  // clang-format on

  inputs_grp->get_option("--input-prefix")->excludes(inputs_grp->get_option("--input-files"));
  params_grp->get_option("--max-memory")
      ->excludes(params_grp->get_option("--max-memory-per-thread"));

  _config = std::monostate{};
}

void Cli::make_metadata_subcommand() {
  [[maybe_unused]] auto &sc =
      *_cli.add_subcommand("metadata",
                           "Fetch the metadata from one of the .parquet file generated by NCHG.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = MetadataConfig{};
           });

  _config = MetadataConfig{};
  auto &c = std::get<MetadataConfig>(_config);

  // clang-format off
  sc.add_option(
    "parquet",
    c.input_path,
    "Path a .parquet file generated by NCHG.")
    ->check(IsValidParquet)
    ->required();
  sc.add_flag(
    "--raw",
    c.raw,
    "Print the metadata as is without any decoding.")
    ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

void Cli::make_view_subcommand() {
  [[maybe_unused]] auto &sc =
      *_cli.add_subcommand("view", "View records stored in a .parquet file.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = ViewConfig{};
           });

  _config = ViewConfig{};
  auto &c = std::get<ViewConfig>(_config);

  // clang-format off
  sc.add_option(
    "parquet",
    c.input_path,
    "Path to the .parquet file to be viewed.")
    ->check(IsValidParquet)
    ->required();
  sc.add_option(
    "--range",
    c.range1,
    "Coordinates of the genomic regions to be fetched following UCSC-style notation (chr1:0-1000).")
    ->capture_default_str();
  sc.add_option(
    "--range2",
    c.range2,
    "Coordinates of the genomic regions to be fetched following UCSC-style notation (chr1:0-1000).")
    ->capture_default_str();
  sc.add_option(
    "--pvalue-cutoff",
    c.pvalue_cutoff,
    "P-value cutoff used to filter records.")
    ->check(CLI::Bound(0.0, 1.0))
    ->capture_default_str();
  sc.add_option(
    "--log-ratio-cutoff",
    c.log_ratio_cutoff,
    "Log-ratio cutoff used to filter records.")
    ->capture_default_str();
  sc.add_flag(
    "--write-header,!--no-write-header",
    c.with_header,
    "Write the file header to stdout.")
    ->capture_default_str();

  // clang-format on

  _config = std::monostate{};
}

void Cli::validate_args() const {
  using enum subcommand;
  switch (get_subcommand()) {
    case cartesian_product:
      validate_cartesian_product_subcommand();  // NOLINT
      break;
    case checksum:
      validate_checksum_subcommand();  // NOLINT
      break;
    case compute:
      validate_compute_subcommand();  // NOLINT
      break;
    case expected:
      validate_expected_subcommand();  // NOLINT
      break;
    case filter:
      validate_filter_subcommand();  // NOLINT
      break;
    case merge:
      validate_merge_subcommand();  // NOLINT
      break;
    case metadata:
      validate_metadata_subcommand();  // NOLINT
      break;
    case view:
      validate_view_subcommand();  // NOLINT
      break;
    case none:
      break;
  }
}

void Cli::validate_cartesian_product_subcommand() const {}

void Cli::validate_checksum_subcommand() const {}

void Cli::validate_compute_subcommand() const {
  const auto &c = std::get<ComputePvalConfig>(_config);
  const auto &sc = *_cli.get_subcommand("compute");

  std::vector<std::string> errors;

  const auto is_mcool = hictk::cooler::utils::is_multires_file(c.path_to_hic.string());
  const auto is_hic = hictk::hic::utils::is_hic_file(c.path_to_hic.string());
  const auto resolution_parsed = !sc.get_option("--resolution")->empty();

  if ((is_mcool || is_hic) && !resolution_parsed) {
    errors.emplace_back(
        "--resolution is a mandatory argument when the input file is in .hic or .mcool format.");
  }

  if (c.output_prefix.parent_path().empty() && !c.chrom1.has_value()) {
    assert(!c.chrom2.has_value());
    errors.emplace_back(
        fmt::format("invalid output prefix \"{}\": the output prefix should have a folder "
                    "component (e.g. \"output_folder/prefix\")",
                    c.output_prefix));
  }

  if (c.min_delta >= c.max_delta) {
    errors.emplace_back("--min-delta should be less than --max-delta.");
  }

  const auto min_delta_parsed = !sc.get_option("--min-delta")->empty();
  const auto max_delta_parsed = !sc.get_option("--max-delta")->empty();
  const auto trans_only_parsed = !sc.get_option("--trans-only")->empty();

  if (trans_only_parsed && (min_delta_parsed || max_delta_parsed)) {
    _warnings.emplace_back("--min-delta and --max-delta are ignored when --trans-only=true");
  }

  if (!c.path_to_domains.empty() && !sc.get_option("--interaction-aggregation-strategy")->empty()) {
    _warnings.emplace_back(
        "--interaction-aggregation-strategy is ignored when no domains have been specified through "
        "the --domains option");
  }

  if (c.compression_method == "lz4" && c.compression_lvl > 9) {
    _warnings.emplace_back("compression method lz4 supports compression levels up to 9");
  }

  if (c.threads > 1 && c.chrom1.has_value()) {
    _warnings.emplace_back(
        "number of threads set with --threads is ignored because --chrom1 has been specified: "
        "concurrency will be limited to a single thread");
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("the following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_expected_subcommand() const {
  const auto &c = std::get<ExpectedConfig>(_config);

  std::vector<std::string> errors;

  if (!c.force && std::filesystem::exists(c.output_path)) {
    errors.emplace_back(
        fmt::format("Refusing to overwrite file {}. Pass --force to overwrite.", c.output_path));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("the following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_filter_subcommand() const {
  const auto &c = std::get<FilterConfig>(_config);

  std::vector<std::string> errors;
  if (!c.force && std::filesystem::exists(c.output_path)) {
    errors.emplace_back(
        fmt::format("Refusing to overwrite file {}. Pass --force to overwrite.", c.output_path));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("the following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_merge_subcommand() const {
  const auto &c = std::get<MergeConfig>(_config);
  const auto &sc = *_cli.get_subcommand("merge");

  std::vector<std::string> errors;

  const auto min_memory = IsValidMemoryDuckDB.lb() * c.threads;
  if (c.memory_limit.value_or(min_memory) < min_memory) {
    // NOLINTBEGIN(*-unchecked-optional-access)
    errors.emplace_back(fmt::format(
        "--max-memory should be at least {} per thread, found {} ({} per thread)",
        pretty_format_memory(IsValidMemoryDuckDB.lb()), pretty_format_memory(*c.memory_limit),
        pretty_format_memory(*c.memory_limit / c.threads)));
    // NOLINTEND(*-unchecked-optional-access)
  }

  if (!sc.get_option("--input-files")->empty() && !sc.get_option("--ignore-report-file")->empty()) {
    _warnings.emplace_back(
        "--ignore-report-file is ignored when files are specified through the --input-files "
        "option");
  }

  if (c.compression_method == "lz4" && c.compression_lvl > 9) {
    _warnings.emplace_back("compression method lz4 supports compression levels up to 9");
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("the following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_metadata_subcommand() const {}

void Cli::validate_view_subcommand() const {}

void Cli::transform_args() {
  using enum subcommand;
  switch (get_subcommand()) {
    case cartesian_product:
      transform_args_cartesian_product_subcommand();  // NOLINT
      break;
    case checksum:
      transform_args_checksum_subcommand();  // NOLINT
      break;
    case compute:
      transform_args_compute_subcommand();  // NOLINT
      break;
    case expected:
      transform_args_expected_subcommand();  // NOLINT
      break;
    case filter:
      transform_args_filter_subcommand();  // NOLINT
      break;
    case merge:
      transform_args_merge_subcommand();  // NOLINT
      break;
    case metadata:
      transform_args_metadata_subcommand();  // NOLINT
      break;
    case view:
      transform_args_view_subcommand();  // NOLINT
      break;
    case none:
      break;
  }
}

void Cli::transform_args_cartesian_product_subcommand() {
  auto &c = std::get<CartesianProductConfig>(_config);

  if (c.chrom1.has_value()) {
    if (!c.chrom2.has_value()) {
      c.chrom2 = c.chrom1;
    }
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity <= SPDLOG_LEVEL_CRITICAL);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

void Cli::transform_args_checksum_subcommand() {
  auto &c = std::get<ChecksumConfig>(_config);

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity <= SPDLOG_LEVEL_CRITICAL);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

static std::string get_path_to_executable() {
#ifdef _WIN32
  std::string path(MAX_PATH, '\0');
  if (GetModuleFileNameA(NULL, path.data(), path.size())) {
    return std::string{path.c_str()};
  }

#elif defined(__APPLE__)
  std::string path(PATH_MAX, '\0');
  std::uint32_t count = PATH_MAX;
  if (!_NSGetExecutablePath(path.data(), &count)) {
    return path.substr(0, count);
  }

#else
  std::string path(PATH_MAX, '\0');
  const auto count = readlink("/proc/self/exe", path.data(), path.size());
  if (count != -1) {
    return path.substr(0, static_cast<std::size_t>(count));
  }
#endif
  throw std::runtime_error("unable to generate the path to NCHG.");
}

void Cli::transform_args_compute_subcommand() {
  auto &c = std::get<ComputePvalConfig>(_config);
  if (c.chrom1.has_value()) {
    if (!c.chrom2.has_value()) {
      c.chrom2 = c.chrom1;
    }

    if (c.chrom1 == c.chrom2) {
      c.compute_cis = true;
      c.compute_trans = false;
    } else {
      c.compute_cis = false;
      c.compute_trans = true;
    }

    c.output_path = c.output_prefix;
    c.output_prefix.clear();

    if (c.output_path.extension() != ".parquet") {
      c.output_path = std::filesystem::path{fmt::format("{}.parquet", c.output_path.string())};
    }

    c.threads = 1;
  }

  c.exec = get_path_to_executable();

  if (c.compression_method == "lz4") {
    c.compression_lvl = std::min(c.compression_lvl, std::uint8_t{9});
  }

  // NOLINTNEXTLINE(*-mt-unsafe)
  if (const auto *path = std::getenv("NCHG_LOG_MESSAGE_QUEUE_NAME"); path) {
    c.log_message_queue = path;
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity <= SPDLOG_LEVEL_CRITICAL);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

void Cli::transform_args_expected_subcommand() {
  auto &c = std::get<ExpectedConfig>(_config);
  if (c.chrom1 != "all" && c.chrom2 == "all") {
    c.chrom2 = c.chrom1;
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity <= SPDLOG_LEVEL_CRITICAL);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

void Cli::transform_args_filter_subcommand() {
  auto &c = std::get<FilterConfig>(_config);
  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity <= SPDLOG_LEVEL_CRITICAL);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

void Cli::transform_args_merge_subcommand() {
  auto &c = std::get<MergeConfig>(_config);

  if (!c.memory_limit.has_value()) {
    c.memory_limit = c.memory_per_thread * c.threads;
  }

  if (c.compression_method == "lz4") {
    c.compression_lvl = std::min(c.compression_lvl, std::uint8_t{9});
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity <= SPDLOG_LEVEL_CRITICAL);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

void Cli::transform_args_metadata_subcommand() {}

void Cli::transform_args_view_subcommand() {}

bool Cli::handle_help_flags() {
  if (_help_flag.empty()) {
    return false;
  }

  if (_help_flag == "license") {
    fmt::print("{}", config::license::license);
    /*
    } else if (_help_flag == "cite") {
      fmt::print("{}", get_citation());
    */
  } else if (_help_flag == "docs") {
    fmt::println("https://github.com/paulsengroup/NCHG?tab=readme-ov-file#NCHG");
  } else {
    unreachable_code();
  }

  _subcommand = subcommand::none;
  _exit_code = 0;
  return true;
}

}  // namespace nchg
