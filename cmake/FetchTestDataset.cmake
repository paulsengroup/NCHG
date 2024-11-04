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
  DOWNLOAD https://zenodo.org/records/14036359/files/nchg_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=a556c5d0022167d10247a8657f8d0299f3609cd3c11b78071d91ebd28ae753f3
  "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst"
)
# gersemi: on

file(ARCHIVE_EXTRACT INPUT "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst" DESTINATION "${PROJECT_SOURCE_DIR}")
