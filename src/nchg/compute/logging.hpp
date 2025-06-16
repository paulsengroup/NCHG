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

#pragma once

#include <spdlog/spdlog.h>

#include <atomic>
#include <boost/interprocess/ipc/message_queue.hpp>
#include <cstddef>
#include <filesystem>
#include <glaze/beve.hpp>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>

namespace glz {

template <>
struct from<BEVE, spdlog::log_clock::time_point> {
  template <auto Opts>
  static void op(spdlog::log_clock::time_point &timestamp, auto &&...args) {
    using BuffT = decltype(timestamp.time_since_epoch().count());

    BuffT buff{};
    parse<BEVE>::op<Opts>(buff, args...);
    timestamp = spdlog::log_clock::time_point(spdlog::log_clock::duration{buff});
  }
};

template <>
struct from<BEVE, spdlog::string_view_t> {
  template <auto Opts>
  static void op(spdlog::string_view_t &payload, auto &&...args) {
    std::string_view buff{};
    parse<BEVE>::op<Opts>(buff, args...);
    payload = spdlog::string_view_t{buff.data(), buff.size()};
  }
};

template <>
struct to<BEVE, spdlog::log_clock::time_point> {
  template <auto Opts>
  static void op(const spdlog::log_clock::time_point &timestamp, auto &&...args) {
    serialize<BEVE>::op<Opts>(timestamp.time_since_epoch().count(), args...);
  }
};

template <>
struct to<BEVE, spdlog::string_view_t> {
  template <auto Opts>
  static void op(spdlog::string_view_t s, auto &&...args) {
    serialize<BEVE>::op<Opts>(std::string_view{s.data(), s.size()}, args...);
  }
};

template <>
struct meta<spdlog::details::log_msg> {
  using T = spdlog::details::log_msg;
  static constexpr auto value = object("t", &T::time, "l", &T::level, "p", &T::payload);
};

}  // namespace glz

namespace nchg {

class MessageQueue {
  std::string _name;
  std::unique_ptr<boost::interprocess::message_queue> _queue;
  std::mutex _mtx;
  std::string _eoq_signal;
  const std::atomic<bool> *_early_return{};
  bool _destroy_queue{false};

  static constexpr std::size_t _max_message_size = sizeof(char) * 2048UL;

  MessageQueue(std::string name, std::size_t num_proc, const std::atomic<bool> *early_return,
               bool create);

 public:
  MessageQueue() = default;

  [[nodiscard]] static MessageQueue create(const std::string &name, std::size_t num_proc,
                                           const std::atomic<bool> &early_return);
  [[nodiscard]] static MessageQueue open(const std::string &name);

  MessageQueue(const MessageQueue &other) = delete;
  MessageQueue(MessageQueue &&other) noexcept;

  ~MessageQueue() noexcept;

  MessageQueue &operator=(const MessageQueue &other) = delete;
  MessageQueue &operator=(MessageQueue &&other) noexcept;

  [[nodiscard]] std::string_view name() const noexcept;

  void send(const spdlog::details::log_msg &msg);
  void send(std::string_view msg);
  void send(std::string_view msg, std::unique_lock<std::mutex> lck);

  bool receive();
  void send_eoq_signal();
  void close() noexcept;

 private:
  [[nodiscard]] static std::unique_ptr<boost::interprocess::message_queue> open_queue(
      const std::filesystem::path &path);
  [[nodiscard]] static std::unique_ptr<boost::interprocess::message_queue> create_queue(
      const std::filesystem::path &path, std::size_t num_proc);
};

bool setup_file_backed_logger(const std::filesystem::path &output_prefix, bool force);

void setup_beve_logger(const std::string &name);

}  // namespace nchg
