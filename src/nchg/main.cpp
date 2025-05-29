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

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <deque>
#include <memory>
#include <mutex>
#include <tuple>
#include <utility>

#include "nchg/common.hpp"
#include "nchg/tools/cli.hpp"
#include "nchg/tools/config.hpp"
#include "nchg/tools/tools.hpp"
#include "nchg/version.hpp"

using namespace nchg;

template <std::size_t CAPACITY>
class GlobalLogger {
  //                                   [2021-08-12 17:49:34.581] [info]: my log msg
  static constexpr auto *_msg_pattern{"[%Y-%m-%d %T.%e] %^[%l]%$: %v"};
  using NCHGLogMsg = std::pair<spdlog::level::level_enum, std::string>;
  std::deque<NCHGLogMsg> _msg_buffer;

  std::mutex _mtx;
  std::atomic<std::size_t> _num_msg_enqueued{0};
  std::atomic<bool> _ok{false};

  [[nodiscard]] static std::shared_ptr<spdlog::sinks::stderr_color_sink_mt> init_stderr_sink() {
    constexpr spdlog::level::level_enum level{SPDLOG_ACTIVE_LEVEL};
    auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    stderr_sink->set_pattern(_msg_pattern);
    stderr_sink->set_level(level);

    return stderr_sink;
  }

  [[nodiscard]] std::shared_ptr<spdlog::sinks::callback_sink_mt> init_warning_callback_sink() {
    if constexpr (CAPACITY != 0 && SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_WARN) {
      auto callback_sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
          [this](const spdlog::details::log_msg &msg) noexcept { enqueue_msg(msg); });
      callback_sink->set_pattern(_msg_pattern);
      callback_sink->set_level(spdlog::level::warn);

      return callback_sink;
    }

    return {};
  }

  template <typename... T>
  void print_noexcept(fmt::format_string<T...> fmt, T &&...args) noexcept {
    try {
      fmt::print(stderr, fmt, std::forward<T>(args)...);
    } catch (...) {  // NOLINT
      _ok = false;
    }
  }

  void enqueue_msg(const spdlog::details::log_msg &msg) noexcept {
    if (msg.level < spdlog::level::warn) [[likely]] {
      return;
    }

    ++_num_msg_enqueued;

    try {
      [[maybe_unused]] const std::scoped_lock lck(_mtx);
      if (_msg_buffer.size() == CAPACITY) [[unlikely]] {
        _msg_buffer.pop_front();
      }
      _msg_buffer.emplace_back(msg.level, std::string{msg.payload.begin(), msg.payload.end()});
    } catch (...) {  // NOLINT
    }
  }

  void replay_warnings() {
    [[maybe_unused]] const std::scoped_lock lck(_mtx);
    if (_msg_buffer.empty()) {
      return;
    }
    auto stderr_sink = init_stderr_sink();
    stderr_sink->set_level(spdlog::level::warn);

    auto logger = std::make_shared<spdlog::logger>("tmp_logger", std::move(stderr_sink));
    logger->set_level(spdlog::level::warn);

    if (_num_msg_enqueued <= _msg_buffer.size()) {
      logger->warn("replaying the last {} warning message(s)", _num_msg_enqueued.load());
    } else {
      logger->warn("replaying the last {}/{} warning messages", _msg_buffer.size(),
                   _num_msg_enqueued.load());
    }
    for (const auto &msg : _msg_buffer) {
      logger->log(msg.first, msg.second);
    }
    _msg_buffer.clear();
  }

  static void reset_logger() noexcept {
    try {
      spdlog::set_default_logger(
          std::make_shared<spdlog::logger>("main_logger", init_stderr_sink()));
    } catch (...) {  // NOLINT
    }
  }

 public:
  GlobalLogger() noexcept {
    try {
      spdlog::set_default_logger(std::make_shared<spdlog::logger>(
          "main_logger",
          spdlog::sinks_init_list{init_stderr_sink(), init_warning_callback_sink()}));
      _ok = true;
    } catch (const std::exception &e) {
      print_noexcept("FAILURE! Failed to setup NCHG's logger: {}", e.what());
    } catch (...) {
      print_noexcept("FAILURE! Failed to setup NCHG's logger: unknown error");
    }
  }

  GlobalLogger(const GlobalLogger &other) = delete;
  GlobalLogger(GlobalLogger &&other) noexcept = default;

  GlobalLogger &operator=(const GlobalLogger &other) = delete;
  GlobalLogger &operator=(GlobalLogger &&other) noexcept = default;

  ~GlobalLogger() noexcept {
    if (!_ok) {
      reset_logger();
      return;
    }

    try {
      replay_warnings();
    } catch (const std::exception &e) {
      print_noexcept("FAILURE! Failed to replay NCHG warnings: {}", e.what());
    } catch (...) {
      print_noexcept("FAILURE! Failed to replay NCHG warnings: unknown error");
    }
    reset_logger();
  }

  static void set_level(int lvl) {
    if (auto logger = spdlog::default_logger(); logger) {
      for (auto &sink : logger->sinks()) {
        sink->set_level(std::max(sink->level(), spdlog::level::level_enum{lvl}));
      }
      logger->set_level(spdlog::level::level_enum{lvl});
    }
  }

  void print_welcome_msg() {
    if (_ok) {
      SPDLOG_INFO("Running {}", config::version::str_long());
    }
  }

  [[nodiscard]] constexpr bool ok() const noexcept { return _ok; }

  void clear() noexcept {
    if (_ok) {
      _msg_buffer.clear();
      _num_msg_enqueued = 0;
    }
  }
};

// NOLINTNEXTLINE(*-err58-cpp, *-avoid-non-const-global-variables)
static auto global_logger = std::make_unique<GlobalLogger<256>>();

static auto acquire_global_logger() noexcept { return std::move(global_logger); }

[[nodiscard]] static bool should_print_welcome_msg(Cli::subcommand subcmd, const Config &config) {
  if (subcmd == Cli::subcommand::none) {
    return false;
  }

  if (std::holds_alternative<ComputePvalConfig>(config)) {
    return std::get<ComputePvalConfig>(config).log_message_queue.empty();
  }

  return true;
}

[[nodiscard]] static bool should_print_cli_warnings(const Config &config) {
  const auto *c = std::get_if<ComputePvalConfig>(&config);

  if (!c) {
    return true;
  }

  if (!c->log_message_queue.empty()) {
    return false;
  }

  return true;
}

template <typename Logger>
static std::tuple<int, Cli::subcommand, Config> parse_cli_and_setup_logger(Cli &cli,
                                                                           Logger &logger) {
  try {
    auto config = cli.parse_arguments();
    const auto subcmd = cli.get_subcommand();
    const auto ec = cli.exit();
    std::visit(
        [&]<typename T>(const T &config_) {
          if constexpr (!std::is_same_v<T, std::monostate>) {
            if (logger.ok()) {
              logger.set_level(config_.verbosity);  // NOLINT
              if (should_print_welcome_msg(subcmd, config)) {
                logger.print_welcome_msg();
              }
            }
          }
        },
        config);

    return std::make_tuple(ec, subcmd, config);
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli.exit(e), Cli::subcommand::none, Config{});
  } catch (const std::filesystem::filesystem_error &e) {
    SPDLOG_ERROR("FAILURE! {}", e.what());
    return std::make_tuple(1, Cli::subcommand::none, Config());
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(stderr,
               "FAILURE! An error occurred while setting up the main "
               "application logger: {}.\n",
               e.what());
    return std::make_tuple(1, Cli::subcommand::none, Config());
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  std::unique_ptr<Cli> cli{nullptr};

  try {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    auto local_logger = acquire_global_logger();
    cli = std::make_unique<Cli>(argc, argv);
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(*cli, *local_logger);
    if (ec != 0 || subcmd == Cli::subcommand::none) {
      local_logger->clear();
      return ec;
    }

    if (should_print_cli_warnings(config)) {
      cli->log_warnings();
    }

    assert(!config.valueless_by_exception());
    return std::visit(
        []<typename Config>(const Config &c) -> int {
          if constexpr (std::is_same_v<Config, std::monostate>) {
            throw std::runtime_error(
                "Caught attempt to visit a Config object in std::monostate. "
                "This should never be possible! "
                "If you see this message, please file an issue on GitHub.");
          } else {
            return run_command(c);
          }
        },
        config);

  } catch (const CLI::ParseError &e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error
                          //  messages (if any)
  } catch (const std::bad_alloc &err) {
    SPDLOG_CRITICAL("FAILURE! Unable to allocate enough memory: {}", err.what());
    return 1;
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(stderr, "FAILURE! NCHG encountered the following error while logging: {}\n",
               e.what());
    return 1;
  } catch (const std::exception &e) {
    if (cli) {
      SPDLOG_CRITICAL("FAILURE! NCHG {} encountered the following error: {}",
                      cli->get_printable_subcommand(), e.what());
    } else {
      SPDLOG_CRITICAL("FAILURE! NCHG encountered the following error: {}", e.what());
    }
    return 1;
  } catch (...) {
    SPDLOG_CRITICAL(
        "FAILURE! NCHG encountered the following error: caught an unhandled exception!\n"
        "If you see this message, please file an issue on GitHub.");
    return 1;
  }
  return 0;
}
