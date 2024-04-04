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

#include <fmt/format.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdio>
#include <tuple>

#include "nchg/cli.hpp"
#include "nchg/config.hpp"
#include "nchg/tools.hpp"

using namespace nchg;

static void setup_logger_console() {
  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  //                        [2021-08-12 17:49:34.581] [info]: my log msg
  stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");

  auto main_logger = std::make_shared<spdlog::logger>("main_logger", stderr_sink);

  spdlog::set_default_logger(main_logger);
}

static void setup_logger_console(int verbosity_lvl, bool print_version) {
  for (auto &sink : spdlog::default_logger()->sinks()) {
    sink->set_level(spdlog::level::level_enum(verbosity_lvl));
  }

  if (print_version) {
    SPDLOG_INFO(FMT_STRING("Running nchg v{}"), "0.0.1");
  }
}

static std::tuple<int, Cli::subcommand, Config> parse_cli_and_setup_logger(Cli &cli) {
  try {
    auto config = cli.parse_arguments();
    const auto subcmd = cli.get_subcommand();
    std::visit(
        [&](const auto &config_) {
          using T = hictk::remove_cvref_t<decltype(config_)>;
          constexpr auto is_monostate = std::is_same_v<T, std::monostate>;
          if constexpr (!is_monostate) {
            setup_logger_console(config_.verbosity, subcmd != Cli::subcommand::help);
          }
        },
        config);

    return std::make_tuple(0, subcmd, config);
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli.exit(e), Cli::subcommand::help, Config());
  } catch (const std::filesystem::filesystem_error &e) {
    SPDLOG_ERROR(FMT_STRING("FAILURE! {}"), e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! An error occurred while setting up the main "
                          "application logger: {}.\n"),
               e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  std::unique_ptr<Cli> cli{nullptr};
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  try {
    setup_logger_console();
    cli = std::make_unique<Cli>(argc, argv);
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(*cli);
    if (ec != 0 || subcmd == Cli::subcommand::help) {
      return ec;
    }

    using sc = Cli::subcommand;
    switch (subcmd) {
      case sc::compute: {
        return run_nchg_compute(std::get<ComputePvalConfig>(config));
      }
      case sc::expected: {
        return run_nchg_expected(std::get<ExpectedConfig>(config));
      }
      case sc::filter: {
        return run_nchg_filter(std::get<FilterConfig>(config));
      }
      case sc::help:  // NOLINT
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in nchg::main() should be unreachable! "
            "If you see this message, please file an issue on GitHub");
    }

  } catch (const CLI::ParseError &e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error
                          //  messages (if any)
  } catch (const std::bad_alloc &err) {
    fmt::print(stderr, FMT_STRING("FAILURE! Unable to allocate enough memory: {}\n"), err.what());
    return 1;
  } catch (const std::exception &e) {
    if (cli) {
      fmt::print(stderr, FMT_STRING("FAILURE! nchg encountered the following error: {}\n"),
                 e.what());
    } else {
      fmt::print(stderr, FMT_STRING("FAILURE! hictk encountered the following error: {}\n"),
                 e.what());
    }
    return 1;
  } catch (...) {
    fmt::print(stderr, FMT_STRING("FAILURE! nchg encountered the following error: Caught an "
                                  "unhandled exception! "
                                  "If you see this message, please file an issue on GitHub.\n"));
    return 1;
  }
  return 0;
}
