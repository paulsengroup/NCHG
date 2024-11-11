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
#include <hictk/transformers/join_genomic_coords.hpp>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "nchg/concepts.hpp"
#include "nchg/expected_matrix.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {
namespace internal {

template <typename N>
  requires arithmetic<N>
[[nodiscard]] constexpr double compute_odds_ratio(N n, N total_sum, N sum1, N sum2) {
  if constexpr (std::is_floating_point_v<N>) {
    if (std::isnan(n)) [[unlikely]] {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  if (sum1 == 0 || sum2 == 0) [[unlikely]] {
    return std::numeric_limits<double>::quiet_NaN();
  }

  assert(n <= sum1);
  assert(n <= sum2);
  assert(sum1 <= total_sum);
  assert(sum2 <= total_sum);
  assert(sum1 + sum2 <= total_sum + n);

  const auto num = n * (total_sum - sum1 - sum2 + n);
  const auto denom = (sum1 - n) * (sum2 - n);

  assert(num >= 0);
  assert(denom >= 0);

  return conditional_static_cast<double>(num) / conditional_static_cast<double>(denom);
}

// GCC gets confused if we declare this as a member function
// (static, not static, inline, constexpr, noexcept(true), noexcept(false)... It does not matter)
template <typename N1,
          typename N2 = std::conditional_t<std::is_floating_point_v<N1>, double, std::uint64_t>>
  requires arithmetic<N1>
constexpr std::pair<double, N2> aggregate_marginals(const hictk::GenomicInterval &range,
                                                    std::uint32_t resolution,
                                                    const std::vector<N1> &marginals,
                                                    const std::vector<bool> &bin_mask) noexcept {
  assert(resolution > 0);

  const auto i1 = range.start() / resolution;
  const auto i2 = (range.end() + resolution - 1) / resolution;

  if (i1 == i2) [[unlikely]] {
    return {N2{}, N2{}};
  }

  assert(i2 <= marginals.size());
  assert(i2 <= bin_mask.size());

  std::size_t bin_masked = 0;
  N2 sum{};
  for (auto i = i1; i < i2; ++i) {
    sum += conditional_static_cast<N2>(marginals[i]);
    bin_masked += bin_mask[i];
  }

  const auto bin_masked_frac = static_cast<double>(bin_masked) / static_cast<double>(i2 - i1);

  return {bin_masked_frac, sum};
}

[[nodiscard]] inline std::pair<std::uint64_t, double> aggregate_pixels_cis(
    const hictk::File &f, const ExpectedMatrixStats &expected_matrix,
    const hictk::GenomicInterval &range1, const hictk::GenomicInterval &range2,
    std::uint64_t min_delta, std::uint64_t max_delta, const std::vector<bool> &bin_mask) {
  assert(range1.chrom() == range2.chrom());
  assert(min_delta <= max_delta);

  std::uint64_t obs{};
  double exp{};

  const auto query_offset1 =
      static_cast<std::int64_t>(f.bins().at(range1.chrom(), range1.start()).rel_id());
  const auto query_offset2 =
      static_cast<std::int64_t>(f.bins().at(range2.chrom(), range2.start()).rel_id());

  const auto query_height = (range1.size() + f.resolution() - 1) / f.resolution();
  const auto query_width = (range2.size() + f.resolution() - 1) / f.resolution();

  std::visit(
      [&](const auto &f_) {
        const auto sel = f_.fetch(range1.chrom().name(), range1.start(), range1.end(),
                                  range2.chrom().name(), range2.start(), range2.end());

        const hictk::transformers::JoinGenomicCoords jsel(
            sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), f.bins_ptr());

        for (const hictk::Pixel<std::uint32_t> &p : jsel) {
          const auto delta = p.coords.bin2.start() - p.coords.bin1.start();

          const auto i1 = p.coords.bin1.rel_id();
          const auto i2 = p.coords.bin2.rel_id();

          if (delta < min_delta || delta >= max_delta || bin_mask[i1] || bin_mask[i2])
              [[unlikely]] {
            continue;
          }

          const auto observed_count = p.count;
          const auto expected_count = expected_matrix.at(i1, i2);

          const auto j1 = static_cast<std::int64_t>(i1) - query_offset1;
          const auto j2 = static_cast<std::int64_t>(i2) - query_offset2;
          const auto j3 = static_cast<std::int64_t>(i2) - query_offset1;
          const auto j4 = static_cast<std::int64_t>(i1) - query_offset2;

          bool added = false;
          if (j1 >= 0 && j1 < query_height && j2 >= 0 && j2 < query_width) {
            obs += observed_count;
            exp += expected_count;
            added = true;
          }

          if (added && j1 == j3 && j2 == j4) [[unlikely]] {
            continue;
          }

          if (j3 >= 0 && j3 < query_height && j4 >= 0 && j4 < query_width) {
            obs += observed_count;
            exp += expected_count;
          }
        }
      },
      f.get());

  return {obs, exp};
}

[[nodiscard]] inline std::pair<std::uint64_t, double> aggregate_pixels_trans(
    const hictk::File &f, const ExpectedMatrixStats &expected_matrix,
    const hictk::GenomicInterval &range1, const hictk::GenomicInterval &range2,
    const std::vector<bool> &bin_mask1, const std::vector<bool> &bin_mask2) {
  assert(range1.chrom() != range2.chrom());

  std::uint64_t obs{};
  double exp{};

  const auto chrom_offset1 = f.bins().at(range1.chrom(), 0).id();
  const auto chrom_offset2 = f.bins().at(range2.chrom(), 0).id();

  std::visit(
      [&](const auto &f_) {
        const auto sel = f_.fetch(range1.chrom().name(), range1.start(), range1.end(),
                                  range2.chrom().name(), range2.start(), range2.end());

        std::for_each(sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(),
                      [&](const hictk::ThinPixel<std::uint32_t> &p) {
                        assert(p.bin1_id >= chrom_offset1);
                        assert(p.bin2_id >= chrom_offset2);
                        const auto i1 = p.bin1_id - chrom_offset1;
                        const auto i2 = p.bin2_id - chrom_offset2;

                        if (!bin_mask1[i1] && !bin_mask2[i2]) [[likely]] {
                          obs += p.count;
                          exp += expected_matrix.at(i1, i2);
                        }
                      });
      },
      f.get());

  return {obs, exp};
}

}  // namespace internal

inline std::uint64_t NCHG::compute_N1(const hictk::GenomicInterval &range1,
                                      const hictk::GenomicInterval &range2,
                                      double max_bad_bin_threshold) const noexcept {
  assert(_obs_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] = internal::aggregate_marginals(
      range1, _fp->resolution(), _obs_matrix->marginals1(),
      *_expected_values.bin_mask(range1.chrom(), range2.chrom()).first);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<std::uint64_t>::max();
  }

  return sum;
}

inline std::uint64_t NCHG::compute_N2(const hictk::GenomicInterval &range1,
                                      const hictk::GenomicInterval &range2,
                                      double max_bad_bin_threshold) const noexcept {
  assert(_obs_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] = internal::aggregate_marginals(
      range2, _fp->resolution(), _obs_matrix->marginals2(),
      *_expected_values.bin_mask(range1.chrom(), range2.chrom()).second);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<std::uint64_t>::max();
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

  if (range1.chrom() == range2.chrom()) {
    const auto &mask = _expected_values.bin_mask(range1.chrom());
    assert(mask);
    const auto [obs, exp] = internal::aggregate_pixels_cis(
        *_fp, *_exp_matrix, range1, range2, params().min_delta, params().max_delta, *mask);

    return Result{obs, exp};
  }

  const auto &[mask1, mask2] = _expected_values.bin_mask(range1.chrom(), range2.chrom());

  assert(mask1);
  assert(mask2);

  const auto [obs, exp] =
      internal::aggregate_pixels_trans(*_fp, *_exp_matrix, range1, range2, *mask1, *mask2);

  return Result{obs, exp};
}

template <typename N>
  requires arithmetic<N>
inline NCHGResult NCHG::compute_stats(hictk::Pixel<N> pixel, double exp, N obs_sum, N N1, N N2,
                                      double exp_sum, double L1, double L2,
                                      std::vector<double> &buffer) {
  constexpr auto bad_sum = std::numeric_limits<N>::max();
  if (pixel.count == 0 || N1 == bad_sum || N2 == bad_sum) [[unlikely]] {
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

  if (!intra_matrix) {
    obs_sum *= 2;
  }

  const auto odds_ratio = internal::compute_odds_ratio(pixel.count, obs_sum, N1, N2);
  const auto omega = intra_matrix ? internal::compute_odds_ratio(exp, exp_sum, L1, L2) : 1.0;
  const auto log_ratio = std::log2(static_cast<double>(pixel.count)) - std::log2(exp);

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

  if (!intra_matrix) {
    obs_sum /= 2;
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
