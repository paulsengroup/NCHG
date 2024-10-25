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

include(FetchContent)

set(HICTK_ENABLE_TESTING OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_WITH_EIGEN OFF)
set(HICTK_BUILD_TOOLS OFF)
set(HICTK_INSTALL OFF)
FetchContent_Declare(
        hictk
        GIT_REPOSITORY https://github.com/paulsengroup/hictk.git
        GIT_TAG        v2.0.1
        EXCLUDE_FROM_ALL
        OVERRIDE_FIND_PACKAGE
        SYSTEM)

FetchContent_Declare(
        stocc
        URL "${CMAKE_CURRENT_SOURCE_DIR}/../../external/stocc-20240322.tar.xz"
        URL_HASH SHA256=91e414f44445dcdfa569148c2484a50cadd786e30d833836d98f6111f82e9850
        UPDATE_DISCONNECTED 1
        EXCLUDE_FROM_ALL
        SYSTEM)

find_package(FMT REQUIRED)
find_package(hictk REQUIRED)
find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)


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
  INTERFACE fmt::fmt-header-only
            hictk::libhictk
            phmap
            spdlog::spdlog_header_only
            stocc)

set_target_properties(libnchg PROPERTIES OUTPUT_NAME libnchg)
