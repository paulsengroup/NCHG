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

find_package(Arrow REQUIRED)
find_package(bshoshany-thread-pool REQUIRED)
find_package(FMT REQUIRED)
find_package(glaze REQUIRED)
find_package(hictk REQUIRED)
find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)
find_package(Threads REQUIRED)
find_package(xxHash REQUIRED)

add_library(nchg_io STATIC)
add_library(nchg::io ALIAS nchg_io)

target_sources(
  nchg_io
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/file_hashing.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_metadata.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_store.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/k_merger_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/parquet_stats_file_helpers.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/parquet_stats_file_reader.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/parquet_stats_file_reader_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/parquet_stats_file_writer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/parquet_stats_file_writer_impl.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/text.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/text_impl.hpp"
)

target_include_directories(
  nchg_io
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  nchg_io
  PRIVATE
    nchg_project_options
    nchg_project_warnings
    nchg::version
    fmt::fmt-header-only
    spdlog::spdlog_header_only
  PUBLIC
    nchg::common
    Arrow::arrow_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
    bshoshany-thread-pool::bshoshany-thread-pool
    glaze::glaze
    hictk::libhictk
    Parquet::parquet_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
    phmap
    Threads::Threads
    xxHash::xxhash
)
