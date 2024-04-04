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

#pragma once

#include "nchg/config.hpp"

namespace nchg {

[[nodiscard]] int run_nchg_compute(const ComputePvalConfig& c);
[[nodiscard]] int run_nchg_filter(const FilterConfig& c);
[[nodiscard]] int run_nchg_expected(const ExpectedConfig& c);

}  // namespace nchg
