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
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(HighFive REQUIRED)
find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)

add_library(nchg_descriptive_stats STATIC)
add_library(nchg::descriptive_stats ALIAS nchg_descriptive_stats)

target_sources(
  nchg_descriptive_stats
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/expected_matrix.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/expected_matrix_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/expected_values.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/expected_values_aggregator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/expected_values_aggregator_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/expected_values_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix_stats_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/observed_matrix.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/observed_matrix_impl.hpp"
)

target_include_directories(
  nchg_descriptive_stats
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  nchg_descriptive_stats
  PRIVATE
    nchg_project_options
    nchg_project_warnings
    fmt::fmt-header-only
    spdlog::spdlog_header_only
  PUBLIC
    nchg::common
    nchg::preproc
    hictk::libhictk
    HDF5::HDF5
    HighFive
    phmap
)
