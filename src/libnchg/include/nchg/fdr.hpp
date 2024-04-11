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

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include "nchg/common.hpp"

namespace nchg {

template <typename Stats>
class BH_FDR {
  std::vector<Stats> _pvalues{};
  std::vector<std::size_t> _idx{};
  std::vector<std::size_t> _ranks{};

 public:
  explicit BH_FDR(std::vector<Stats> pvalues_);

  template <typename UnaryOperation = identity>
  [[nodiscard]] auto correct(UnaryOperation op = identity()) -> std::vector<Stats>;
};

}  // namespace nchg

#include "../../fdr_impl.hpp"
