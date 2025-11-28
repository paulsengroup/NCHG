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

find_package(Boost REQUIRED)
find_package(FMT REQUIRED)
find_package(hictk REQUIRED)

add_library(genomic_domains STATIC)
add_library(nchg::genomic_domains ALIAS genomic_domains)

target_sources(
  genomic_domains
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/genomic_domains.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/genomic_domains_impl.hpp"
)

target_include_directories(genomic_domains PUBLIC include/)
target_link_libraries(
  genomic_domains
  PRIVATE
    nchg_project_options
    nchg_project_warnings
  PUBLIC
    Boost::headers
    nchg::common
    hictk::libhictk
)
