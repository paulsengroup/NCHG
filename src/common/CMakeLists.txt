# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: GPL-3.0
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public License as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public License along
# with this library.  If not, see
# <https://www.gnu.org/licenses/>.

find_package(FMT REQUIRED)
find_package(hictk REQUIRED)
find_package(spdlog REQUIRED)

add_library(nchg_common INTERFACE)
add_library(nchg::common ALIAS nchg_common)

target_sources(
  nchg_common
  INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include/nchg/common.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/nchg/concepts.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/nchg/hash.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/nchg/suppress_warnings.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/nchg/type_traits.hpp"
)

target_include_directories(
  nchg_common
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  nchg_common
  INTERFACE
    nchg_project_options
    nchg_project_warnings
    fmt::fmt-header-only
    hictk::libhictk
    spdlog::spdlog_header_only
)
