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
// <https://www.gnu.org/licenses/>.

#include "./logging.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/spdlog.h>

#include <atomic>
#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/interprocess/timed_utils.hpp>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <glaze/beve.hpp>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

#include "nchg/version.hpp"

namespace nchg {

MessageQueue::MessageQueue(std::string name, std::size_t num_proc,
                           const std::atomic<bool> *early_return, bool create)
    : _name(std::move(name)),
      _queue(create ? create_queue(_name, num_proc) : open_queue(_name)),
      _eoq_signal(
          fmt::format("### EOQ {}", std::chrono::system_clock::now().time_since_epoch().count())),
      _early_return(early_return),
      _destroy_queue(create) {}

MessageQueue MessageQueue::create(const std::string &name, std::size_t num_proc,
                                  const std::atomic<bool> &early_return) {
  return {name, num_proc, &early_return, true};
}

MessageQueue MessageQueue::open(const std::string &name) { return {name, 1, nullptr, false}; }

MessageQueue::MessageQueue(MessageQueue &&other) noexcept
    : _name(std::move(other._name)),
      _queue(std::move(other._queue)),
      _eoq_signal(std::move(other._eoq_signal)),
      _early_return(other._early_return) {}

MessageQueue::~MessageQueue() noexcept { close(); }

MessageQueue &MessageQueue::operator=(MessageQueue &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  _name = std::move(other._name);
  _queue = std::move(other._queue);
  _eoq_signal = std::move(other._eoq_signal);
  _early_return = other._early_return;

  return *this;
}

std::string_view MessageQueue::name() const noexcept { return _name; }

void MessageQueue::send(const spdlog::details::log_msg &msg) {
  assert(_queue);
  std::string buff{};
  std::unique_lock lck(_mtx);
  if (!glz::write_beve(msg, buff)) [[likely]] {
    send(buff, std::move(lck));
  }
}

void MessageQueue::send(std::string_view msg) { send(msg, std::unique_lock{_mtx}); }

void MessageQueue::send(std::string_view msg, std::unique_lock<std::mutex> lck) {
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

bool MessageQueue::receive() {
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
        msg.time, msg.source, msg.level != spdlog::level::info ? msg.level : spdlog::level::trace,
        msg.payload);
  }

  return true;
}

void MessageQueue::send_eoq_signal() { send(_eoq_signal); }

void MessageQueue::close() noexcept {
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

std::unique_ptr<boost::interprocess::message_queue> MessageQueue::open_queue(
    const std::filesystem::path &path) {
  SPDLOG_DEBUG("opening message queue with name={}...", path);
  return std::make_unique<boost::interprocess::message_queue>(boost::interprocess::open_only,
                                                              path.c_str());
}

std::unique_ptr<boost::interprocess::message_queue> MessageQueue::create_queue(
    const std::filesystem::path &path, std::size_t num_proc) {
  assert(num_proc != 0);
  SPDLOG_DEBUG("initializing a message queue with name={}...", path);
  boost::interprocess::message_queue::remove(path.c_str());
  return std::make_unique<boost::interprocess::message_queue>(
      boost::interprocess::create_only, path.c_str(), num_proc * 4, sizeof(char) * 2048UL);
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

bool setup_file_backed_logger(const std::filesystem::path &output_prefix, bool force) {
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
    file_sink->set_level(spdlog::level::trace);
    file_sink->set_pattern(log_pattern);

    auto tmp_logger = std::make_shared<spdlog::logger>("tmp_logger", file_sink);
    tmp_logger->set_level(spdlog::level::info);
    tmp_logger->log(spdlog::level::info,
                    fmt::format("log file was generated by {}", config::version::str_long()));

    logger->sinks().emplace_back(std::move(file_sink));
    logger->set_level(spdlog::level::trace);

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

void setup_beve_logger(const std::string &name) {
  auto logger = spdlog::default_logger();
  if (!logger || name.empty()) {
    return;
  }

  auto message_queue = std::make_shared<MessageQueue>(MessageQueue::open(name));

  auto callback_sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
      [message_queue](const spdlog::details::log_msg &msg) { message_queue->send(msg); });

  callback_sink->set_level(spdlog::level::trace);

  logger->sinks().clear();
  logger->sinks().emplace_back(std::move(callback_sink));
  logger->set_level(spdlog::level::trace);
}

}  // namespace nchg
