# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
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

function(ConfigureLicensing license_file input_config_file output_config_file)
  set(PRE_CONFIGURE_FILE "${input_config_file}")
  set(POST_CONFIGURE_FILE "${output_config_file}")

  file(READ "${license_file}" NCHG_LICENSE)

  if(NOT WIN32)
    file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
  endif()

  configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
endfunction()
