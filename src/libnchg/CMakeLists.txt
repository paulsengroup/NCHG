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
  URL_HASH SHA256=06ca370cfa93b3619d379adc8e016f28a0dce0dd644fbce6d28118b8ce889059
  SYSTEM)

set(NCHG_BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS OFF)
FetchContent_MakeAvailable(stocc)
set(BUILD_SHARED_LIBS "${NCHG_BUILD_SHARED_LIBS}")
unset(NCHG_BUILD_SHARED_LIBS)

add_library(libnchg INTERFACE)
add_library(nchg::libnchg ALIAS libnchg)

target_sources(
  libnchg
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(libnchg INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                             "$<INSTALL_INTERFACE:include>")

target_link_libraries(
  libnchg
  INTERFACE Boost::headers
            fmt::fmt-header-only
            hictk::libhictk
            phmap
            spdlog::spdlog_header_only
            stocc)

set_target_properties(libnchg PROPERTIES OUTPUT_NAME libnchg)
