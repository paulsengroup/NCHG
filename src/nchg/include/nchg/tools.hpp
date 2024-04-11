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

#ifndef _WIN32
#include <csignal>
#endif

#include <atomic>
#include <cstddef>

#include "nchg/config.hpp"

namespace nchg {

#ifdef _WIN32
[[nodiscard]] int run_nchg_compute(const ComputePvalConfig& c, std::atomic<std::uint32_t*>& pids,
                                   const std::atomic<std::size_t>& num_pids);
#else
[[nodiscard]] int run_nchg_compute(const ComputePvalConfig& c, std::atomic<pid_t*>& pids,
                                   const std::atomic<std::size_t>& num_pids);
#endif

[[nodiscard]] int run_nchg_filter(const FilterConfig& c);
[[nodiscard]] int run_nchg_expected(const ExpectedConfig& c);

}  // namespace nchg
