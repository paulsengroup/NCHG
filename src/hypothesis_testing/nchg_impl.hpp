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

template <typename PixelSelector>
inline NCHG::iterator<PixelSelector>::iterator(PixelSelector selector,
                                               std::shared_ptr<const ObservedMatrix> obs,
                                               std::shared_ptr<const ExpectedMatrixStats> exp,
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
      _max_delta(max_delta),
      _read_value(true) {
  jump_to_next_valid_pixel();

  if (_pixel_it != _sentinel_it) {
    std::ignore = **this;
  }
}

template <typename PixelSelector>
inline auto NCHG::iterator<PixelSelector>::at_end(PixelSelector selector,
                                                  std::shared_ptr<const ObservedMatrix> obs,
                                                  std::shared_ptr<const ExpectedMatrixStats> exp)
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
      _value((other._value)),
      _read_value(other._read_value) {}

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
  assert(_pixel_it != _sentinel_it);
  assert(_sel);
  assert(_obs);
  assert(_exp);
  assert(_buffer);

  if (!_read_value) {
    return _value;
  }

  const auto intra_matrix = _obs->chrom1() == _obs->chrom2();

  const auto obs_sum = _obs->sum();
  const auto exp_sum = _exp->sum();

  // clang-format off
  hictk::Pixel p{
    _sel->bins().at(_pixel_it->bin1_id),
    _sel->bins().at(_pixel_it->bin2_id),
    static_cast<std::uint64_t>(_pixel_it->count)
  };
  // clang-format on

  const auto N1 = compute_N1(p);
  const auto N2 = compute_N2(p);
  const auto L1 = compute_L1(p);
  const auto L2 = compute_L2(p);

  const auto exp = _exp->at(p.coords.bin1.rel_id(), p.coords.bin2.rel_id());

  const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
  if (intra_matrix && (delta < _min_delta || delta >= _max_delta)) [[unlikely]] {
    _read_value = false;
    _value = {std::move(p), exp, 1.0, 0.0, 0.0, 0.0};
    return _value;
  }

  _value = compute_stats(std::move(p), exp, obs_sum, N1, N2, exp_sum, L1, L2, *_buffer);
  _read_value = false;

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
  _read_value = true;
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

    assert(bin1_id < _bin_mask1->size());
    assert(bin2_id < _bin_mask2->size());

    const auto bin1_masked = (*_bin_mask1)[bin1_id];
    const auto bin2_masked = (*_bin_mask2)[bin2_id];

    if (!bin1_masked && !bin2_masked) [[likely]] {
      break;
    }
    std::ignore = ++_pixel_it;
  }
  _read_value = true;
}

template <typename PixelSelector>
constexpr std::uint64_t NCHG::iterator<PixelSelector>::compute_N1(
    const hictk::Pixel<std::uint64_t> &pixel) const noexcept {
  assert(_obs);
  const auto i = pixel.coords.bin1.rel_id();
  assert(i < _obs->marginals1().size());
  return _obs->marginals1()[i];
}

template <typename PixelSelector>
constexpr std::uint64_t NCHG::iterator<PixelSelector>::compute_N2(
    const hictk::Pixel<std::uint64_t> &pixel) const noexcept {
  assert(_obs);
  const auto i = pixel.coords.bin2.rel_id();
  assert(i < _obs->marginals2().size());
  return _obs->marginals2()[i];
}

template <typename PixelSelector>
constexpr double NCHG::iterator<PixelSelector>::compute_L1(
    const hictk::Pixel<std::uint64_t> &pixel) const noexcept {
  assert(_exp);
  const auto i = pixel.coords.bin1.rel_id();
  assert(i < _exp->marginals1().size());
  return _exp->marginals1()[i];
}

template <typename PixelSelector>
constexpr double NCHG::iterator<PixelSelector>::compute_L2(
    const hictk::Pixel<std::uint64_t> &pixel) const noexcept {
  assert(_exp);
  const auto i = pixel.coords.bin2.rel_id();
  assert(i < _exp->marginals2().size());
  return _exp->marginals2()[i];
}

}  // namespace nchg
