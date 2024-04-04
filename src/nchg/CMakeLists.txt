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

find_package(CLI11 REQUIRED)
find_package(FMT REQUIRED)
find_package(hictk REQUIRED)
find_package(HighFive REQUIRED)
find_package(spdlog REQUIRED)

add_executable(NCHG)

target_sources(
  NCHG
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/compute/compute.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/expected/expected.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/filter/filter.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp"
          "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

target_include_directories(NCHG PRIVATE include/)
target_link_libraries(
  NCHG
  PUBLIC libnchg
         CLI11::CLI11
         fmt::fmt-header-only
         hictk::libhictk
         HighFive
         spdlog::spdlog_header_only)
