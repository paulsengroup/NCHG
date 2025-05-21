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
# <https://www.gnu.org/licenses/>

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git" AND NCHG_ENABLE_GIT_VERSION_TRACKING)
  message(
    WARNING
    "-- Unable to find .git/ under \"${PROJECT_SOURCE_DIR}\". Setting -DNCHG_ENABLE_GIT_VERSION_TRACKING=OFF"
  )
  set(NCHG_ENABLE_GIT_VERSION_TRACKING OFF)
endif()

function(ConfigureVersioning input_config_folder output_config_folder)
  set(PRE_CONFIGURE_FILE "${input_config_folder}/git.hpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/git.hpp")

  if(NCHG_ENABLE_GIT_VERSION_TRACKING)
    include(FetchContent)
    FetchContent_Declare(
      _nchg_cmake-git-version-tracking
      URL
        "${PROJECT_SOURCE_DIR}/external/cmake-git-version-tracking.20250308.tar.xz"
      URL_HASH "SHA256=1e3f448fa746dde298d0175a51e83a0149614debdd6f4379e2db0f67acc8d20b"
    )
    FetchContent_MakeAvailable(_nchg_cmake-git-version-tracking)

    set(GIT_IGNORE_UNTRACKED ON)
    include("${_nchg_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake")
  else()
    # Add dummy target
    add_custom_target(_nchg_check_git)

    if(NOT DEFINED NCHG_GIT_RETRIEVED_STATE)
      set(NCHG_GIT_RETRIEVED_STATE false)
    endif()
    if(NOT DEFINED NCHG_GIT_HEAD_SHA1)
      set(NCHG_GIT_HEAD_SHA1 "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_IS_DIRTY)
      set(NCHG_GIT_IS_DIRTY false)
    endif()
    if(NOT DEFINED NCHG_GIT_AUTHOR_NAME)
      set(NCHG_GIT_AUTHOR_NAME "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_AUTHOR_EMAIL)
      set(NCHG_GIT_AUTHOR_EMAIL "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_COMMIT_DATE_ISO8601)
      set(NCHG_GIT_COMMIT_DATE_ISO8601 "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_COMMIT_SUBJECT)
      set(NCHG_GIT_COMMIT_SUBJECT "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_COMMIT_BODY)
      set(NCHG_GIT_COMMIT_BODY "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_DESCRIBE)
      set(NCHG_GIT_DESCRIBE "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_BRANCH)
      set(NCHG_GIT_BRANCH "unknown")
    endif()
    if(NOT DEFINED NCHG_GIT_TAG)
      set(NCHG_GIT_TAG "unknown")
    endif()

    if(NOT WIN32)
      file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
    endif()
    configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
  endif()

  set(PRE_CONFIGURE_FILE "${input_config_folder}/version.hpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/version.hpp")

  if(NOT WIN32)
    file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
  endif()
  configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
endfunction()
