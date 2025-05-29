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
# <https://www.gnu.org/licenses/>

if(NOT WIN32)
  file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)
endif()

set(TEST_DATASET_TAR "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst")

message(STATUS "Fetching the test dataset")

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/15546722/files/nchg_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=5337b8c117492b3f3bec68aef8406e89b81ed6741ebaf42b41bb43d5e161b3b6
  "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst"
)
# gersemi: on

message(STATUS "Fetching the test dataset - done")

message(STATUS "Extracting the test dataset")

file(ARCHIVE_EXTRACT INPUT "${TEST_DATASET_TAR}" DESTINATION "${PROJECT_SOURCE_DIR}")

message(STATUS "Extracting the test dataset - done")
message(STATUS "Test datasets can be found under \"${PROJECT_SOURCE_DIR}/test/data/\"")
