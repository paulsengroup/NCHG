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
set(HICTK_ENABLE_GIT_VERSION_TRACKING OFF)

FetchContent_Declare(
  hictk
  GIT_REPOSITORY https://github.com/paulsengroup/hictk.git
  GIT_TAG v2.1.1
  EXCLUDE_FROM_ALL
  OVERRIDE_FIND_PACKAGE
  SYSTEM
)

add_subdirectory(common)
add_subdirectory(descriptive_stats)
add_subdirectory(hypothesis_testing)
add_subdirectory(io)
add_subdirectory(preproc)

add_subdirectory(nchg)
