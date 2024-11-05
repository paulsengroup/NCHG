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
#include <fmt/std.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdio>
#include <deque>
#include <mutex>
#include <tuple>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/tools/cli.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"

using namespace nchg;

static constexpr auto *nchg_log_msg_pattern{"[%Y-%m-%d %T.%e] %^[%l]%$: %v"};
static std::mutex nchg_global_mtx;
static constexpr std::size_t nchg_warning_message_buffer_capacity = 256;
using NCHGLogMsg = std::pair<spdlog::level::level_enum, std::string>;
static std::deque<NCHGLogMsg> nchg_warning_message_buffer{};

static void setup_logger_console() {
  //                   [2021-08-12 17:49:34.581] [info]: my log msg

  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  auto callback_sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [](const spdlog::details::log_msg &msg) mutable {
        if (msg.level == spdlog::level::warn) [[unlikely]] {
          [[maybe_unused]] const std::scoped_lock lck(nchg_global_mtx);
          if (nchg_warning_message_buffer.size() == nchg_warning_message_buffer_capacity)
              [[unlikely]] {
            nchg_warning_message_buffer.pop_front();
          }
          nchg_warning_message_buffer.emplace_back(
              msg.level, std::string{msg.payload.begin(), msg.payload.end()});
        }
      });

  stderr_sink->set_pattern(nchg_log_msg_pattern);
  callback_sink->set_pattern(nchg_log_msg_pattern);
  stderr_sink->set_level(spdlog::level::debug);
  callback_sink->set_level(spdlog::level::warn);

  auto main_logger = std::make_shared<spdlog::logger>(
      "main_logger", spdlog::sinks_init_list{std::move(stderr_sink), std::move(callback_sink)});

  spdlog::set_default_logger(std::move(main_logger));
}

static void setup_logger_console(int verbosity_lvl, bool print_version) {
  for (auto &sink : spdlog::default_logger()->sinks()) {
    sink->set_level(std::max(sink->level(), spdlog::level::level_enum{verbosity_lvl}));
  }
  spdlog::set_level(spdlog::level::level_enum{verbosity_lvl});

  if (print_version) {
    SPDLOG_INFO("Running NCHG v{}", "0.0.2");
  }
}

static std::tuple<int, Cli::subcommand, Config> parse_cli_and_setup_logger(Cli &cli) {
  try {
    auto config = cli.parse_arguments();
    const auto subcmd = cli.get_subcommand();
    const auto ec = cli.exit();
    std::visit(
        [&]<typename T>(const T &config_) {
          constexpr auto is_monostate = std::is_same_v<T, std::monostate>;
          if constexpr (!is_monostate) {
            setup_logger_console(config_.verbosity, subcmd != Cli::subcommand::help);
          }
        },
        config);

    return std::make_tuple(ec, subcmd, config);
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli.exit(e), Cli::subcommand::help, Config());
  } catch (const std::filesystem::filesystem_error &e) {
    SPDLOG_ERROR("FAILURE! {}", e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(stderr,
               "FAILURE! An error occurred while setting up the main "
               "application logger: {}.\n",
               e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  }
}

template <typename Callable, typename Config>
[[nodiscard]] static int run_subcommand(Callable fx, const Config &config) {
  auto replay_warnings = [] {
    [[maybe_unused]] const std::scoped_lock lck(nchg_global_mtx);
    if (nchg_warning_message_buffer.empty()) {
      return;
    }
    auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_st>();
    stderr_sink->set_pattern(nchg_log_msg_pattern);
    stderr_sink->set_level(spdlog::level::warn);

    auto logger = std::make_shared<spdlog::logger>("tmp_logger", std::move(stderr_sink));
    logger->set_level(spdlog::level::warn);

    logger->warn("replaying the last {} warning message(s)", nchg_warning_message_buffer.size());
    for (const auto &msg : nchg_warning_message_buffer) {
      logger->log(msg.first, msg.second);
    }
    nchg_warning_message_buffer.clear();
  };

  try {
    const auto ec = fx(config);
    replay_warnings();
    return ec;
  } catch (...) {
    replay_warnings();
    throw;
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

    cli->log_warnings();

    using sc = Cli::subcommand;
    switch (subcmd) {
      case sc::compute: {
        return run_subcommand(run_nchg_compute, std::get<ComputePvalConfig>(config));
      }
      case sc::expected: {
        return run_subcommand(run_nchg_expected, std::get<ExpectedConfig>(config));
      }
      case sc::filter: {
        return run_subcommand(run_nchg_filter, std::get<FilterConfig>(config));
      }
      case sc::merge: {
        return run_subcommand(run_nchg_merge, std::get<MergeConfig>(config));
      }
      case sc::view: {
        return run_subcommand(run_nchg_view, std::get<ViewConfig>(config));
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
    SPDLOG_CRITICAL("FAILURE! Unable to allocate enough memory: {}\n", err.what());
    return 1;
  } catch (const std::exception &e) {
    if (cli) {
      SPDLOG_CRITICAL("FAILURE! NCHG {} encountered the following error: {}\n",
                      cli->get_printable_subcommand(), e.what());
    } else {
      SPDLOG_CRITICAL("FAILURE! NCHG encountered the following error: {}\n", e.what());
    }
    return 1;
  } catch (...) {
    SPDLOG_CRITICAL(
        "FAILURE! NCHG encountered the following error: caught an "
        "unhandled exception!\n"
        "If you see this message, please file an issue on GitHub.\n");
    return 1;
  }
  return 0;
}
