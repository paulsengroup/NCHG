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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <hictk/pixel.hpp>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {
namespace internal {
[[nodiscard]] constexpr double compute_odds_ratio(double n, double total_sum, double sum1,
                                                  double sum2) {
  if (std::isnan(n) || sum1 == 0 || sum2 == 0) [[unlikely]] {
    return std::numeric_limits<double>::quiet_NaN();
  }

  assert(n <= sum1);
  assert(n <= sum2);
  assert(sum1 <= total_sum);
  assert(sum2 <= total_sum);
  assert(sum1 + sum2 <= total_sum);

  const auto num = (n * (total_sum - sum1 - sum2 + n));
  const auto denom = (sum1 - n) * (sum2 - n);

  return num / denom;
}
}  // namespace internal

template <typename PixelSelector>
inline NCHG::iterator<PixelSelector>::iterator(PixelSelector selector,
                                               std::shared_ptr<const ObservedMatrix> obs,
                                               std::shared_ptr<const ExpectedMatrix> exp,
                                               std::shared_ptr<const std::vector<bool>> bin_mask1,
                                               std::shared_ptr<const std::vector<bool>> bin_mask2,
                                               std::uint64_t min_delta, std::uint64_t max_delta)
    : _sel(std::make_shared<const PixelSelector>(std::move(selector))),
      _pixel_it(_sel->template begin<N>()),
      _sentinel_it(_sel->template end<N>()),
      _obs(std::move(obs)),
      _exp(std::move(exp)),
      _bin_mask1(std::move(bin_mask1)),
      _bin_mask2(std::move(bin_mask2)),
      _min_delta(min_delta),
      _max_delta(max_delta) {
  jump_to_next_valid_pixel();
}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::at_end(PixelSelector selector,
                                                  std::shared_ptr<const ObservedMatrix> obs,
                                                  std::shared_ptr<const ExpectedMatrix> exp)
    -> iterator {
  iterator it{};
  it._sel = std::make_shared<const PixelSelector>(std::move(selector));
  it._pixel_it = it._sel->template end<N>();
  it._sentinel_it = it._pixel_it;
  it._obs = std::move(obs);
  it._exp = std::move(exp);

  return it;
}

template <typename PixelSelector>
inline NCHG::iterator<PixelSelector>::iterator(const iterator &other)
    : _sel(other._sel),
      _pixel_it(other._pixel_it),
      _sentinel_it(other._sentinel_it),
      _obs(other._obs),
      _exp(other._exp),
      _bin_mask1(other._bin_mask1),
      _bin_mask2(other._bin_mask2),
      _min_delta(other._min_delta),
      _max_delta(other._max_delta),
      _value((other._value)) {}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::operator=(const iterator &other) -> iterator & {
  if (this == &other) {
    return *this;
  }

  _sel = other._sel;
  _pixel_it = other._pixel_it;
  _sentinel_it = other._sentinel_it;
  _obs = other._obs;
  _exp = other._exp;
  _bin_mask1 = other._bin_mask1;
  _bin_mask2 = other._bin_mask2;
  _min_delta = other._min_delta;
  _max_delta = other._max_delta;
  _value = other._value;

  return *this;
}

template <typename PixelSelector>
inline bool NCHG::iterator<PixelSelector>::operator==(const iterator &other) const noexcept {
  return _pixel_it == other._pixel_it && _obs == other._obs && _exp == other._exp;
}

template <typename PixelSelector>
inline bool NCHG::iterator<PixelSelector>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename PixelSelector>
inline bool NCHG::iterator<PixelSelector>::operator<(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it < other._pixel_it;
}

template <typename PixelSelector>
inline bool NCHG::iterator<PixelSelector>::operator<=(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it <= other._pixel_it;
}

template <typename PixelSelector>
inline bool NCHG::iterator<PixelSelector>::operator>(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it > other._pixel_it;
}

template <typename PixelSelector>
inline bool NCHG::iterator<PixelSelector>::operator>=(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it >= other._pixel_it;
}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::operator*() const -> const_reference {
  const auto &obs_marginals1 = _obs->marginals1();
  const auto &obs_marginals2 = _obs->marginals2();
  const auto &exp_marginals1 = _exp->marginals1();
  const auto &exp_marginals2 = _exp->marginals2();

  const auto intra_matrix = _obs->chrom1() == _obs->chrom2();

  const auto obs_sum = _obs->sum();
  const auto exp_sum = _exp->sum();

  constexpr double cutoff = 1.0e-20;

  const hictk::Pixel p{_sel->bins(), *_pixel_it};
  const auto i1 = p.coords.bin1.rel_id();
  const auto i2 = p.coords.bin2.rel_id();

  const auto N1 = static_cast<double>(obs_marginals1[i1]);
  const auto N2 = static_cast<double>(obs_marginals2[i2]);
  const auto obs = static_cast<double>(p.count);

  const auto L1 = exp_marginals1[i1];
  const auto L2 = exp_marginals2[i2];

  const auto exp = std::max(_exp->at(i1, i2), cutoff);

  const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
  if (intra_matrix && (delta < _min_delta || delta >= _max_delta)) [[unlikely]] {
    _value = {p, exp, 1.0, 0.0, 0.0, 0.0};
    return _value;
  }

  const auto odds_ratio = internal::compute_odds_ratio(obs, static_cast<double>(obs_sum), N1, N2);
  const auto omega = intra_matrix ? internal::compute_odds_ratio(exp, exp_sum, L1, L2) : 1;

  const auto log_ratio = std::log2(odds_ratio) - std::log2(omega);

  if ((L1 - exp) * (L2 - exp) <= cutoff) [[unlikely]] {
    _value = {p, exp, 1.0, log_ratio, odds_ratio, omega};
    return _value;
  }

  if (!std::isfinite(omega) || omega > odds_ratio) [[unlikely]] {
    _value = {p, exp, 1.0, log_ratio, odds_ratio, omega};
    return _value;
  }

  const auto pvalue =
      compute_pvalue_nchg(*_buffer, static_cast<std::uint64_t>(obs), static_cast<std::uint64_t>(N1),
                          static_cast<std::uint64_t>(N2), obs_sum, omega);
  _value = {p, exp, pvalue, log_ratio, odds_ratio, omega};

  return _value;
}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::operator->() const -> const_pointer {
  _value = **this;
  return &_value;
}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::operator++() -> iterator & {
  std::ignore = ++_pixel_it;
  jump_to_next_valid_pixel();
  return *this;
}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::operator++(int) -> iterator {
  auto it = *this;
  it._str_buffer = std::make_shared<std::vector<double>>();
  std::ignore = ++(*this);
  return it;
}

template <typename PixelSelector>
inline void NCHG::iterator<PixelSelector>::jump_to_next_valid_pixel() {
  assert(_sel);
  assert(_bin_mask1);
  assert(_bin_mask2);

  while (_pixel_it != _sentinel_it) {
    const hictk::Pixel p{_sel->bins(), *_pixel_it};
    const auto bin1_id = p.coords.bin1.rel_id();
    const auto bin2_id = p.coords.bin2.rel_id();

    const auto bin1_masked = (*_bin_mask1)[bin1_id];
    const auto bin2_masked = (*_bin_mask2)[bin2_id];

    if (!bin1_masked && !bin2_masked) [[likely]] {
      break;
    }
    std::ignore = ++_pixel_it;
  }
}

}  // namespace nchg
