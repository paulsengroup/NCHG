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

// clang-format off
#include "nchg/suppress_warnings.hpp"
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/type_fwd.h>
#include <parquet/platform.h>
NCHG_DISABLE_WARNING_POP
// clang-format on

#include <hictk/reference.hpp>
#include <memory>
#include <string_view>
#include <type_traits>

namespace nchg {

// https://stackoverflow.com/a/16000226
template <typename T, typename = int>
struct has_pval_corrected : std::false_type {};

template <typename T>
struct has_pval_corrected<T, decltype((void)T::pval_corrected, 0)> : std::true_type {};

[[nodiscard]] std::shared_ptr<arrow::Schema> get_schema(const hictk::Reference &chroms);
[[nodiscard]] std::shared_ptr<arrow::Schema> get_schema_padj(const hictk::Reference &chroms);

[[nodiscard]] parquet::Compression::type parse_parquet_compression(std::string_view method);

}  // namespace nchg
