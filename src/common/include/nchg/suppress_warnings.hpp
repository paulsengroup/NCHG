// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public License along
// with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// GCC to MSVC codes: https://github.com/srz-zumix/awesome-cpp-warning

// clang-format off

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

// Defines for MSVC
#ifdef _MSC_VER
    #define NCHG_DISABLE_WARNING_PUSH                      __pragma(warning(push))
    #define NCHG_DISABLE_WARNING_POP                       __pragma(warning(pop))
    #define NCHG_DISABLE_WARNING(warningNumber)            __pragma(warning(disable : warningNumber))

    #define NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS   NCHG_DISABLE_WARNING(4996)
    #define NCHG_DISABLE_WARNING_MACRO_REDEFINED           NCHG_DISABLE_WARNING(4005)
    #define NCHG_DISABLE_MAYBE_UNINITIALIZED
    #define NCHG_DISABLE_WARNING_OLD_STYLE_CAST            NCHG_DISABLE_WARNING(26475)
#endif

// Defines for GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
    #define NCHG_DO_PRAGMA(X)                              _Pragma(#X)
    #define NCHG_DISABLE_WARNING_PUSH                      NCHG_DO_PRAGMA(GCC diagnostic push)
    #define NCHG_DISABLE_WARNING_POP                       NCHG_DO_PRAGMA(GCC diagnostic pop)
    #define NCHG_DISABLE_WARNING(warningName)              NCHG_DO_PRAGMA(GCC diagnostic ignored warningName)

    #define NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS   NCHG_DISABLE_WARNING("-Wdeprecated-declarations")
    #define NCHG_DISABLE_WARNING_OLD_STYLE_CAST            NCHG_DISABLE_WARNING("-Wold-style-cast")
#endif

// Defines for GCC only
#if defined(__GNUC__) && !defined(__clang__)
    #define NCHG_DISABLE_WARNING_MACRO_REDEFINED
    #define NCHG_DISABLE_MAYBE_UNINITIALIZED               NCHG_DISABLE_WARNING("-Wmaybe-uninitialized")
#endif

// Defines for Clang only
#ifdef __clang__
    #define NCHG_DISABLE_WARNING_MACRO_REDEFINED           NCHG_DISABLE_WARNING("-Wmacro-redefined")
    #define NCHG_DISABLE_MAYBE_UNINITIALIZED
#endif

// Defines for unknown/unsupported compilers
#if !defined(_MSC_VER) && !defined(__GNUC__) && !defined(__clang__)
    #define NCHG_DISABLE_WARNING
    #define NCHG_DISABLE_WARNING_PUSH
    #define NCHG_DISABLE_WARNING_POP

    #define NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
    #define NCHG_DISABLE_WARNING_MACRO_REDEFINED
    #define NCHG_DISABLE_MAYBE_UNINITIALIZED
    #define NCHG_DISABLE_WARNING_OLD_STYLE_CAST
#endif

// NOLINTEND(cppcoreguidelines-macro-usage)

// clang-format on
