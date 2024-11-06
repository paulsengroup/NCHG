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
  assert(sum1 + sum2 <= 2 * total_sum);

  const auto num = n * ((2 * total_sum) - sum1 - sum2 + n);
  const auto denom = (sum1 - n) * (sum2 - n);

  assert(num >= 0);
  assert(denom >= 0);

  return num / denom;
}

// GCC gets confused if we declare this as a member function
// (static, not static, inline, constexpr, noexcept(true), noexcept(false)... It does not matter)
template <typename N>
constexpr std::pair<double, double> aggregate_marginals(
    const hictk::GenomicInterval &range, std::uint32_t resolution, const std::vector<N> &marginals,
    const std::vector<bool> &bin_mask) noexcept {
  assert(resolution > 0);

  const auto i1 = range.start() / resolution;
  const auto i2 = (range.end() + resolution - 1) / resolution;

  if (i1 == i2) [[unlikely]] {
    return {0.0, 0.0};
  }

  assert(i2 <= marginals.size());
  assert(i2 <= bin_mask.size());

  std::size_t bin_masked = 0;
  double sum = 0.0;
  for (auto i = i1; i < i2; ++i) {
    sum += static_cast<double>(marginals[i]);
    bin_masked += bin_mask[i];
  }

  const auto bin_masked_frac = static_cast<double>(bin_masked) / static_cast<double>(i2 - i1);

  return {bin_masked_frac, sum};
}
}  // namespace internal

inline double NCHG::compute_N1(const hictk::GenomicInterval &range1,
                               const hictk::GenomicInterval &range2,
                               double max_bad_bin_threshold) const noexcept {
  assert(_obs_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] = internal::aggregate_marginals(
      range1, _fp->resolution(), _obs_matrix->marginals1(),
      *_expected_values.bin_mask(range1.chrom(), range2.chrom()).first);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return sum;
}

inline double NCHG::compute_N2(const hictk::GenomicInterval &range1,
                               const hictk::GenomicInterval &range2,
                               double max_bad_bin_threshold) const noexcept {
  assert(_obs_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] = internal::aggregate_marginals(
      range2, _fp->resolution(), _obs_matrix->marginals2(),
      *_expected_values.bin_mask(range1.chrom(), range2.chrom()).second);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return sum;
}

inline double NCHG::compute_L1(const hictk::GenomicInterval &range1,
                               const hictk::GenomicInterval &range2,
                               double max_bad_bin_threshold) const noexcept {
  assert(_exp_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] = internal::aggregate_marginals(
      range1, _fp->resolution(), _exp_matrix->marginals1(),
      *_expected_values.bin_mask(range1.chrom(), range2.chrom()).first);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return sum;
}

inline double NCHG::compute_L2(const hictk::GenomicInterval &range1,
                               const hictk::GenomicInterval &range2,
                               double max_bad_bin_threshold) const noexcept {
  assert(_exp_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] = internal::aggregate_marginals(
      range2, _fp->resolution(), _exp_matrix->marginals2(),
      *_expected_values.bin_mask(range1.chrom(), range2.chrom()).second);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return sum;
}

inline auto NCHG::aggregate_pixels(const hictk::GenomicInterval &range1,
                                   const hictk::GenomicInterval &range2) const {
  assert(_fp);
  assert(_exp_matrix);

  struct Result {
    std::uint64_t obs{};
    double exp{};
  };

  const auto min_delta = params().min_delta;
  const auto max_delta = params().max_delta;

  const auto intra_matrix = range1.chrom() == range2.chrom();

  const auto &mask1 = *_expected_values.bin_mask(range1.chrom(), range2.chrom()).first;
  const auto &mask2 = *_expected_values.bin_mask(range1.chrom(), range2.chrom()).second;

  std::uint64_t obs{};
  double exp{};

  std::visit(
      [&](const auto &f) {
        const auto sel = f.fetch(range1.chrom().name(), range1.start(), range1.end(),
                                 range2.chrom().name(), range2.start(), range2.end());

        const hictk::transformers::JoinGenomicCoords jsel(
            sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), f.bins_ptr());
        for (const hictk::Pixel<std::uint32_t> &p : jsel) {
          const auto delta =
              intra_matrix ? p.coords.bin2.start() - p.coords.bin1.start() : min_delta;

          const auto i1 = p.coords.bin1.rel_id();
          const auto i2 = p.coords.bin2.rel_id();
          if (delta >= min_delta && delta < max_delta && !mask1[i1] && !mask2[i2]) [[likely]] {
            obs += p.count;
            exp += _exp_matrix->at(i1, i2);
          }
        }
      },
      _fp->get());

  return Result{obs, exp};
}

inline NCHGResult NCHG::compute_stats(hictk::Pixel<std::uint64_t> pixel, double exp,
                                      std::uint64_t obs_sum, double N1, double N2, double exp_sum,
                                      double L1, double L2, std::vector<double> &buffer) {
  // N1 and N2 can be NaN when a pixel/domain is masked
  if (pixel.count == 0 || std::isnan(N1) || std::isnan(N2)) [[unlikely]] {
    return {.pixel = std::move(pixel),
            .expected = exp,
            .pval = 1.0,
            .log_ratio = 0.0,
            .odds_ratio = 0.0,
            .omega = 0.0};
  }

  assert(std::isfinite(exp));
  assert(std::isfinite(L1));
  assert(std::isfinite(L2));
  assert(std::isfinite(exp_sum));

  const auto intra_matrix = pixel.coords.bin1.chrom() == pixel.coords.bin2.chrom();

  const auto obs = static_cast<double>(pixel.count);
  const auto odds_ratio = internal::compute_odds_ratio(obs, static_cast<double>(obs_sum), N1, N2);
  const auto omega = intra_matrix ? internal::compute_odds_ratio(exp, exp_sum, L1, L2) : 1.0;

  const auto log_ratio = std::log2(obs) - std::log2(exp);

  if (!std::isfinite(odds_ratio) || !std::isfinite(omega)) [[unlikely]] {
    return {.pixel = std::move(pixel),
            .expected = exp,
            .pval = 1.0,
            .log_ratio = log_ratio,
            .odds_ratio = odds_ratio,
            .omega = omega};
  }

  if (!std::isfinite(omega) || omega > odds_ratio) [[unlikely]] {
    return {.pixel = std::move(pixel),
            .expected = exp,
            .pval = 1.0,
            .log_ratio = log_ratio,
            .odds_ratio = odds_ratio,
            .omega = omega};
  }

  const auto pvalue = compute_pvalue_nchg(buffer, pixel.count, static_cast<std::uint64_t>(N1),
                                          static_cast<std::uint64_t>(N2), obs_sum, omega);
  return {.pixel = std::move(pixel),
          .expected = exp,
          .pval = pvalue,
          .log_ratio = log_ratio,
          .odds_ratio = odds_ratio,
          .omega = omega};
}

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
constexpr double NCHG::iterator<PixelSelector>::compute_N1(
    const hictk::Pixel<std::uint64_t> &pixel) const noexcept {
  assert(_obs);
  const auto i = pixel.coords.bin1.rel_id();
  assert(i < _obs->marginals1().size());
  return static_cast<double>(_obs->marginals1()[i]);
}

template <typename PixelSelector>
constexpr double NCHG::iterator<PixelSelector>::compute_N2(
    const hictk::Pixel<std::uint64_t> &pixel) const noexcept {
  assert(_obs);
  const auto i = pixel.coords.bin2.rel_id();
  assert(i < _obs->marginals2().size());
  return static_cast<double>(_obs->marginals2()[i]);
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
