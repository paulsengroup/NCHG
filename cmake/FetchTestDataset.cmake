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

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/15462237/files/nchg_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=f90b8ab9b13e4da975042faeac5f625fcc265fc43d1ca4ec38314993be2597bc
  "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst"
)
# gersemi: on

file(ARCHIVE_EXTRACT INPUT "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst" DESTINATION "${PROJECT_SOURCE_DIR}")
