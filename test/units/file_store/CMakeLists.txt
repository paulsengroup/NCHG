# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: GPL-3.0
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public
# License along with this library.  If not, see
# <https://www.gnu.org/licenses/>.

find_package(Boost REQUIRED)
find_package(Catch2 REQUIRED)
find_package(glaze REQUIRED)
include(CTest)
include(Catch)

add_executable(nchg_file_store_tests)

target_sources(
  nchg_file_store_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/file_hashing_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/file_result_metadata_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/file_store_test.cpp
)

target_link_libraries(
  nchg_file_store_tests
  PRIVATE
    nchg_project_warnings
    nchg_project_options
    nchg::version
  PUBLIC
    Boost::headers
    glaze::glaze
    nchg::io
)

target_link_system_libraries(
  nchg_file_store_tests
  PUBLIC
  Catch2::Catch2WithMain
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  nchg_file_store_tests
  WORKING_DIRECTORY
    "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
    --skip-benchmarks
)
