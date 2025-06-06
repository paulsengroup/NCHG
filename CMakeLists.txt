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

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.30)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE OFF CACHE BOOL "Enable 'developer mode'")

if(NOT CMAKE_FIND_PACKAGE_PREFER_CONFIG)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
endif()

set(NCHG_PROJECT_VERSION_MAJOR 0)
set(NCHG_PROJECT_VERSION_MINOR 0)
set(NCHG_PROJECT_VERSION_PATCH 3)
set(NCHG_PROJECT_VERSION_SUFFIX "")

project(
  NCHG
  LANGUAGES
    C
    CXX
  VERSION "${NCHG_PROJECT_VERSION_MAJOR}.${NCHG_PROJECT_VERSION_MINOR}.${NCHG_PROJECT_VERSION_PATCH}"
  HOMEPAGE_URL https://github.com/paulsengroup/NCHG
  # DESCRIPTION "TODO"
)

include(FetchContent)
FetchContent_Declare(
  _nchg_project_options
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.36.6.tar.xz"
  URL_HASH SHA256=b4df2a4107847248d4711ee48477cdf93c1a1c7d51a931937759eb9ffba67af8
  SYSTEM
)
FetchContent_MakeAvailable(_nchg_project_options)

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(
      CMAKE_CONFIGURATION_TYPES
      Debug
      Release
      RelWithDebInfo
      CACHE STRING
      "Enabled build types"
      FORCE
    )
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

include("${_nchg_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")

# dynamic_project_options sets recommended defaults and provides user and developer modes and full GUI support for
# choosing options at configure time

# for more flexibility, look into project_options() macro

# Any default can be overridden set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode set(<feature_name>_USER_DEFAULT
# <value>) - set default for user mode

# Initialize project_options variable related to this project This overwrites `project_options` and sets
# `project_warnings` uncomment the options to enable them:

set(ENABLE_CACHE_DEFAULT ON)
set(ENABLE_COMPILE_COMMANDS_SYMLINK_DEFAULT OFF)
set(ENABLE_CONAN_DEFAULT OFF)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_DEFAULT OFF)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT ON)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

set(ENABLE_SANITIZER_ADDRESS_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_ADDRESS_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_LEAK_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_LEAK_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_POINTER_COMPARE_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_COMPARE_DEVELOPER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_UNDEFINED_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_UNDEFINED_DEVELOPER_DEFAULT ON)

dynamic_project_options(
  PREFIX
  nchg
  # CLANG_WARNINGS # Override the defaults for the CLANG warnings GCC_WARNINGS   # Override the defaults for the GCC
  CPPCHECK_OPTIONS
  --enable=performance,portability,style,warning
  --inline-suppr
  # We cannot act on a bug/missing feature of cppcheck
  --suppress=internalAstError
  # if a file does not have an internalAstError, we get an unmatchedSuppression error
  --suppress=unmatchedSuppression
  --suppress=passedByValue
  --inconclusive
  MSVC_WARNINGS
  "${MSVC_WARNINGS}"
  CLANG_WARNINGS
  "${CLANG_WARNINGS}"
  GCC_WARNINGS
  "${GCC_WARNINGS}"
  CUDA_WARNINGS
  "${CUDA_WARNINGS}"
)

target_compile_features(nchg_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

target_compile_definitions(
  nchg_project_options
  INTERFACE
    # Tweak fmt
    FMT_HEADER_ONLY
    # Tweak spdlog
    "SPDLOG_ACTIVE_LEVEL=$<IF:$<CONFIG:Debug>,SPDLOG_LEVEL_TRACE,SPDLOG_LEVEL_INFO>"
    SPDLOG_CLOCK_COARSE
    SPDLOG_FMT_EXTERNAL
    SPDLOG_NO_THREAD_ID
    SPDLOG_NO_ATOMIC_LEVELS
    # Tweak xxHash
    XXH_INLINE_ALL
)

if(WIN32)
  target_compile_definitions(
    nchg_project_options
    INTERFACE
      NOMINMAX
      _CRT_SECURE_NO_WARNINGS
  )
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(NCHG_ENABLE_GIT_VERSION_TRACKING "Retrieve project version and metadata from git" ON)
option(NCHG_ENABLE_TESTING "Build unit tests" ON)
option(NCHG_USE_PIDFD_OPEN "Use pidfd_open() when available" ON)

add_subdirectory(src)

if(NCHG_ENABLE_TESTING)
  enable_testing()
  message(STATUS "Building unit tests.")
  target_compile_definitions(nchg_project_options INTERFACE NCHG_ENABLE_TESTING)
  add_subdirectory(test)
endif()
