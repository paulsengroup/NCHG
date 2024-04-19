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
//

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace nchg {

template <typename Stats>
inline BH_FDR<Stats>::BH_FDR(std::vector<Stats> pvalues_) : _pvalues(std::move(pvalues_)) {}

template <typename Stats>
inline void BH_FDR<Stats>::add_record(Stats&& s) {
  _pvalues.emplace_back(std::move(s));
}

template <typename Stats>
template <typename StatsIt>
inline void BH_FDR<Stats>::add_records(StatsIt first, StatsIt last) {
  _pvalues.insert(_pvalues.end(), first, last);
}

template <typename Stats>
inline void BH_FDR<Stats>::clear() noexcept {
  _pvalues.clear();
  _idx.clear();
  _ranks.clear();
}

template <typename Stats>
template <typename UnaryOperation>
inline auto BH_FDR<Stats>::correct(UnaryOperation op) -> std::vector<Stats> {
  if (_pvalues.size() < 2) {
    return _pvalues;
  }

  _idx.resize(_pvalues.size());
  std::iota(_idx.begin(), _idx.end(), 0);

  std::sort(_idx.begin(), _idx.end(),
            [&](const auto i1, const auto i2) { return op(_pvalues[i1]) < op(_pvalues[i2]); });

  _ranks.resize(_pvalues.size());
  for (std::size_t i = 0; i < _ranks.size(); ++i) {
    _ranks[_idx[i]] = i;
  }

  std::vector<Stats> adj_pvalues(_pvalues.size());
  for (std::size_t i = 0; i < _ranks.size(); ++i) {
    const auto rank = _ranks[i] + 1;
    const auto pctile_rank = static_cast<double>(_pvalues.size()) / static_cast<double>(rank);
    adj_pvalues[i] = _pvalues[i];
    op(adj_pvalues[i]) = std::clamp(op(_pvalues[i]) * pctile_rank, 0.0, 1.0);
  }

  for (std::size_t i = 1; i < _ranks.size(); ++i) {
    auto& pv0 = adj_pvalues[_idx[i - 1]];
    auto& pv1 = adj_pvalues[_idx[i]];

    if (op(pv1) < op(pv0)) {
      op(pv0) = op(pv1);
    }
  }
  return adj_pvalues;
}

}  // namespace nchg
