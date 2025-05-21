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

#pragma once

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <hictk/cooler/uri.hpp>
#include <hictk/cooler/validation.hpp>
#include <hictk/hic/validation.hpp>
#include <string>
#include <string_view>

#include "nchg/tools/config.hpp"

namespace nchg {

class CoolerFileValidator : public CLI::Validator {
 public:
  inline CoolerFileValidator() : Validator("Cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!hictk::cooler::utils::is_cooler(uri)) {
        if (hictk::cooler::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        if (hictk::cooler::utils::is_scool_file(uri)) {
          return "URI points to a .scool file: " + uri;
        }
        const auto path = hictk::cooler::parse_cooler_uri(uri).file_path;
        if (!std::filesystem::exists(path)) {
          return "No such file: " + path;
        }
        return "Not a valid Cooler: " + uri;
      }
      return "";
    };
  }
};

class MultiresCoolerFileValidator : public CLI::Validator {
 public:
  inline MultiresCoolerFileValidator() : Validator("Multires-cooler") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = hictk::cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!hictk::cooler::utils::is_multires_file(uri)) {
        return "Not a valid multi-resolution cooler: " + uri;
      }
      return "";
    };
  }
};

class SingleCellCoolerFileValidator : public CLI::Validator {
 public:
  inline SingleCellCoolerFileValidator() : Validator("Single-cell-cooler") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = hictk::cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!hictk::cooler::utils::is_scool_file(uri)) {
        return "Not a valid single-cell cooler: " + uri;
      }
      return "";
    };
  }
};

class HiCFileValidator : public CLI::Validator {
 public:
  inline HiCFileValidator() : Validator("HiC") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = hictk::cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!hictk::hic::utils::is_hic_file(path)) {
        return "Not a valid .hic file: " + path;
      }
      return "";
    };
  }
};

// clang-format off
inline const auto IsValidCoolerFile = CoolerFileValidator();                      // NOLINT(cert-err58-cpp)
inline const auto IsValidMultiresCoolerFile = MultiresCoolerFileValidator();      // NOLINT(cert-err58-cpp)
inline const auto IsValidSingleCellCoolerFile = SingleCellCoolerFileValidator();  // NOLINT(cert-err58-cpp)
inline const auto IsValidHiCFile = HiCFileValidator();                            // NOLINT(cert-err58-cpp)
// clang-format on

class Cli {
 public:
  enum class subcommand : std::uint_fast8_t {
    none,
    cartesian_product,
    checksum,
    compute,
    expected,
    filter,
    merge,
    view,
  };

  Cli(int argc, char** argv);
  [[nodiscard]] subcommand get_subcommand() const noexcept;
  [[nodiscard]] std::string_view get_printable_subcommand() const noexcept;
  [[nodiscard]] auto parse_arguments() -> Config;
  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] int exit() const noexcept;
  [[nodiscard]] static std::string_view subcommand_to_str(subcommand s) noexcept;
  void log_warnings() const noexcept;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  Config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::none};
  mutable std::vector<std::string> _warnings{};
  std::string_view _help_flag{};

  void make_cartesian_product_subcommand();
  void make_checksum_subcommand();
  void make_compute_subcommand();
  void make_expected_subcommand();
  void make_filter_subcommand();
  void make_merge_subcommand();
  void make_view_subcommand();
  void make_cli();

  void validate_cartesian_product_subcommand() const;
  void validate_checksum_subcommand() const;
  void validate_compute_subcommand() const;
  void validate_expected_subcommand() const;
  void validate_filter_subcommand() const;
  void validate_merge_subcommand() const;
  void validate_view_subcommand() const;
  void validate_args() const;

  void transform_args_cartesian_product_subcommand();
  void transform_args_checksum_subcommand();
  void transform_args_compute_subcommand();
  void transform_args_expected_subcommand();
  void transform_args_filter_subcommand();
  void transform_args_merge_subcommand();
  void transform_args_view_subcommand();
  void transform_args();

  [[nodiscard]] bool handle_help_flags();
};

}  // namespace nchg
