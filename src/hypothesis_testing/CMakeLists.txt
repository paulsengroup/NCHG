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

FetchContent_Declare(
  stocc
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/../../external/stocc-20240322.tar.xz"
  URL_HASH SHA256=b6bb7f778d77f46035614f5adb88c3a990ba36205886f03ce90f7f3e32ad2a44
  UPDATE_DISCONNECTED 1
  EXCLUDE_FROM_ALL
  SYSTEM
)

find_package(FMT REQUIRED)
find_package(hictk REQUIRED)
find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)

FetchContent_MakeAvailable(stocc)

add_library(nchg_hypothesis_testing STATIC)
add_library(nchg::hypothesis_testing ALIAS nchg_hypothesis_testing)

target_sources(
  nchg_hypothesis_testing
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/fdr_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/nchg.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/nchg_impl.hpp"
)

target_include_directories(
  nchg_hypothesis_testing
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  nchg_hypothesis_testing
  PRIVATE
    nchg_project_options
    nchg_project_warnings
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    stocc
  PUBLIC
    nchg::common
    nchg::descriptive_stats
    nchg::genomic_domains
    nchg::preproc
    hictk::libhictk
    phmap
)

# Disable clang-tidy on randomc and stocc
# This seems to be the only reliable way to do so...
configure_file("${PROJECT_SOURCE_DIR}/cmake/.clang-tidy.in" "${randomc_SOURCE_DIR}/.clang-tidy" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/cmake/.clang-tidy.in" "${stocc_SOURCE_DIR}/.clang-tidy" COPYONLY)
