# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: GPL-3.0
#
# This library is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library.  If not, see
# <https://www.gnu.org/licenses/>.

include(GNUInstallDirs)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/libnchg/include/nchg"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  COMPONENT Libraries)

install(TARGETS stocc randomc COMPONENT Libraries)

install(
  TARGETS libnchg
  EXPORT libnchg-targets
  COMPONENT Libraries
  FILE_SET HEADERS
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  PRIVATE_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  INCLUDES
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(
  EXPORT libnchg-targets
  COMPONENT Libraries
  FILE nchgTargets.cmake
  NAMESPACE nchg::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nchg/")

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/nchgConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/nchgConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nchg/")

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/nchgConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nchg"
  COMPONENT Libraries)

install(
  FILES "${PROJECT_SOURCE_DIR}/LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/nchg/"
  COMPONENT Libraries)
