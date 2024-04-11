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
cmake_policy(VERSION 3.25...3.28)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE
    OFF
    CACHE BOOL "Enable 'developer mode'")

if(NOT CMAKE_FIND_PACKAGE_PREFER_CONFIG)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
endif()

# include(cmake/Versioning.cmake)

project(
  NCHG
  LANGUAGES CXX
  VERSION 0.0.1)

include(FetchContent)
FetchContent_Declare(
  _nchg_project_options
  URL "${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.33.0.tar.xz"
  URL_HASH SHA256=b55dddd6c8af37c35b0bdd90d78088ef05beb423d6a56a55850e33fa3d464675
  SYSTEM)
FetchContent_MakeAvailable(_nchg_project_options)

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(CMAKE_CONFIGURATION_TYPES
        Debug Release RelWithDebInfo
        CACHE STRING "Enabled build types" FORCE)
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES
        ${CMAKE_BUILD_TYPE}
        CACHE STRING "Enabled build types" FORCE)
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
set(ENABLE_DOXYGEN_USER OFF)
set(ENABLE_DOXYGEN_DEVELOPER ON)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT ON)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

set(ENABLE_SANITIZER_ADDRESS_USER OFF)
set(ENABLE_SANITIZER_ADDRESS_DEVELOPER ON)
set(ENABLE_SANITIZER_LEAK_USER OFF)
set(ENABLE_SANITIZER_LEAK_DEVELOPER ON)
set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR_USER OFF)
set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR_DEVELOPER ON)

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
  "${CUDA_WARNINGS}")

target_compile_features(nchg_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

# Tweak fmt
target_compile_definitions(nchg_project_options INTERFACE FMT_HEADER_ONLY FMT_ENFORCE_COMPILE_STRING)
# Tweak spdlog
target_compile_definitions(nchg_project_options INTERFACE SPDLOG_FMT_EXTERNAL)

if(WIN32)
  target_compile_definitions(nchg_project_options INTERFACE NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(NCHG_ENABLE_TESTING "Build unit tests" ON)

add_subdirectory(src)

if(NCHG_ENABLE_TESTING)
  enable_testing()
  message(STATUS "Building unit tests.")
  target_compile_definitions(nchg_project_options INTERFACE NCHG_ENABLE_TESTING)
  add_subdirectory(test)
endif()

include(cmake/Install.cmake)
