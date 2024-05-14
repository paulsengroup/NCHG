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
find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(bshoshany-thread-pool REQUIRED)
find_package(CLI11 REQUIRED)
find_package(concurrentqueue REQUIRED)
find_package(FMT REQUIRED)
find_package(hictk REQUIRED)
find_package(Filesystem REQUIRED)
find_package(HighFive REQUIRED)
find_package(spdlog REQUIRED)
find_package(Threads REQUIRED)

add_executable(NCHG)

target_sources(
  NCHG
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/compute/compute.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/expected/expected.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/filter/filter.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/merge/merge.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

target_include_directories(NCHG PRIVATE include/)
target_link_libraries(
  NCHG
  PRIVATE nchg_project_options
          nchg_project_warnings
          libnchg
          arrow::arrow
          Boost::headers
          Boost::filesystem # Required by boost::process::v2
          bshoshany-thread-pool::bshoshany-thread-pool
          CLI11::CLI11
          concurrentqueue::concurrentqueue
          fmt::fmt-header-only
          hictk::libhictk
          HighFive
          spdlog::spdlog_header_only
          std::filesystem
          Threads::Threads)

target_compile_definitions(NCHG PRIVATE SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)

include(GNUInstallDirs)
install(
  TARGETS NCHG
  COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scratch/
  PRIVATE_HEADER DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scratch/)

install(
  FILES "${PROJECT_SOURCE_DIR}/LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/nchg/"
  COMPONENT Runtime)
