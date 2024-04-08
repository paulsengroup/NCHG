# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: GPL-3.0
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library.  If not, see
# <https://www.gnu.org/licenses/>.

find_package(Boost REQUIRED)
find_package(FMT REQUIRED)
find_package(hictk REQUIRED)
find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)

include(FetchContent)

FetchContent_Declare(
  stocc
  URL "${CMAKE_CURRENT_SOURCE_DIR}/../../external/stocc-20240322.tar.xz"
  URL_HASH SHA256=c35afbd5c34f868086561f8342d723b0570ad412864f0f136a3c93c55f9564a8
  SYSTEM)

FetchContent_MakeAvailable(stocc)

set_property(TARGET stocc PROPERTY COMPILE_WARNING_AS_ERROR OFF)

add_library(libnchg INTERFACE)
add_library(nchg::libnchg ALIAS libnchg)

target_sources(
  libnchg
  INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/expected_matrix_impl.hpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/extrapolating_spline_impl.hpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/nchg_impl.hpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/observed_matrix_impl.hpp")

target_include_directories(libnchg INTERFACE include/)
target_link_libraries(
  libnchg
  INTERFACE Boost::headers
            fmt::fmt-header-only
            hictk::libhictk
            phmap
            spdlog::spdlog_header_only
            stocc)
