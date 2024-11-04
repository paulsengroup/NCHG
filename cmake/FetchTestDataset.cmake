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
  DOWNLOAD https://zenodo.org/records/14014123/files/nchg_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=e0a0cff876fbb49d0daca950f00c2d4985eba94553437d5e3d3147b6cf01f648
  "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst"
)
# gersemi: on

file(ARCHIVE_EXTRACT INPUT "${PROJECT_SOURCE_DIR}/test/data/nchg_test_data.tar.zst" DESTINATION "${PROJECT_SOURCE_DIR}")
