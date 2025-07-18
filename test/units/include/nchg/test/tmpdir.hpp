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

#include <filesystem>

#include "nchg/tmpdir.hpp"

namespace nchg::test {

// NOLINTBEGIN(cert-err58-cpp)
inline const TmpDir testdir{true};
inline const std::filesystem::path datadir{"test/data"};
// NOLINTEND(cert-err58-cpp)

}  // namespace nchg::test
