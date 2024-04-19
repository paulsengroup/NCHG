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

#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <stocc.h>

#include <cassert>
#include <hictk/file.hpp>
#include <hictk/fmt/genomic_interval.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <hictk/transformers/pixel_merger.hpp>

#include "nchg/expected_matrix.hpp"
#include "nchg/mad_max_filter.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {

template <typename File>
inline NCHG<File>::NCHG(std::shared_ptr<const File> f, const hictk::Chromosome &chrom1,
                        const hictk::Chromosome &chrom2, double mad_max_, std::uint64_t min_delta,
                        std::uint64_t max_delta)
    : NCHG(f, chrom1, chrom2,
           ExpectedValues<File>::chromosome_pair(f, chrom1, chrom2, mad_max_, min_delta,
                                                 max_delta)) {}

template <typename File>
inline NCHG<File>::NCHG(std::shared_ptr<const File> f, const hictk::Chromosome &chrom1,
                        const hictk::Chromosome &chrom2,
                        ExpectedValues<File> expected_values) noexcept
    : _fp(std::move(f)),
      _chrom1(chrom1),
      _chrom2(chrom2),
      _exp_matrix(init_exp_matrix(chrom1, chrom2, *_fp, expected_values)),
      _obs_matrix(init_obs_matrix(
          chrom1, chrom2, *_fp, *expected_values.bin_mask(chrom1, chrom2).first,
          *expected_values.bin_mask(chrom1, chrom2).second, expected_values.mad_max(),
          expected_values.min_delta(), expected_values.max_delta())),
      _expected_values(std::move(expected_values)) {}

template <typename File>
inline double NCHG<File>::mad_max() const noexcept {
  return _expected_values.mad_max();
}
template <typename File>
inline std::uint64_t NCHG<File>::min_delta() const noexcept {
  return _expected_values.min_delta();
}
template <typename File>
inline std::uint64_t NCHG<File>::max_delta() const noexcept {
  return _expected_values.max_delta();
}

template <typename File>
inline auto NCHG<File>::observed_matrix() const noexcept -> const ObservedMatrix<PixelIt> & {
  return *_obs_matrix;
}

template <typename File>
inline auto NCHG<File>::expected_matrix() const noexcept -> const ExpectedMatrix<PixelIt> & {
  return *_exp_matrix;
}

[[nodiscard]] static double compute_odds_ratio(double n, double total_sum, double sum1,
                                               double sum2) {
  if (std::isnan(n) || sum1 == 0 || sum2 == 0) {
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

template <typename File>
inline NCHG<File>::iterator::iterator(PixelIt pixel_it, PixelIt sentinel_it,
                                      std::shared_ptr<const ObservedMatrix<PixelIt>> obs,
                                      std::shared_ptr<const ExpectedMatrix<PixelIt>> exp,
                                      std::shared_ptr<const std::vector<bool>> bin_mask1,
                                      std::shared_ptr<const std::vector<bool>> bin_mask2,
                                      std::uint64_t min_delta, std::uint64_t max_delta) noexcept
    : _pixel_it(std::move(pixel_it)),
      _sentinel_it(std::move(sentinel_it)),
      _obs(std::move(obs)),
      _exp(std::move(exp)),
      _bin_mask1(std::move(bin_mask1)),
      _bin_mask2(std::move(bin_mask2)),
      _min_delta(min_delta),
      _max_delta(max_delta) {}

template <typename File>
inline NCHG<File>::iterator::iterator(const iterator &other)
    : _pixel_it(other._pixel_it),
      _sentinel_it(other._sentinel_it),
      _obs(other._obs),
      _exp(other._exp),
      _bin_mask1(other._bin_mask1),
      _bin_mask2(other._bin_mask2),
      _min_delta(other._min_delta),
      _max_delta(other._max_delta),
      _value((other._value)) {}

template <typename File>
inline auto NCHG<File>::iterator::operator=(const iterator &other) -> iterator & {
  if (this == &other) {
    return *this;
  }

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

template <typename File>
inline bool NCHG<File>::iterator::operator==(const iterator &other) const noexcept {
  return _pixel_it == other._pixel_it && _obs == other._obs && _exp == other._exp;
}

template <typename File>
inline bool NCHG<File>::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename File>
inline bool NCHG<File>::iterator::operator<(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it < other._pixel_it;
}

template <typename File>
inline bool NCHG<File>::iterator::operator<=(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it <= other._pixel_it;
}

template <typename File>
inline bool NCHG<File>::iterator::operator>(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it > other._pixel_it;
}

template <typename File>
inline bool NCHG<File>::iterator::operator>=(const iterator &other) const noexcept {
  assert(_obs == other._obs);
  assert(_exp == other._exp);

  return _pixel_it >= other._pixel_it;
}

template <typename File>
inline auto NCHG<File>::iterator::operator*() const -> const_reference {
  const auto &obs_marginals1 = _obs->marginals1();
  const auto &obs_marginals2 = _obs->marginals2();
  const auto &exp_marginals1 = _exp->marginals1();
  const auto &exp_marginals2 = _exp->marginals2();

  const auto intra_matrix = _obs->chrom1() == _obs->chrom2();

  const auto obs_sum = _obs->sum();
  const auto exp_sum = _exp->sum();

  const double cutoff = 1.0e-20;

  const auto &p = *_pixel_it;
  const auto i1 = p.coords.bin1.rel_id();
  const auto i2 = p.coords.bin2.rel_id();

  const auto N1 = static_cast<double>(obs_marginals1[i1]);
  const auto N2 = static_cast<double>(obs_marginals2[i2]);
  const auto obs = static_cast<double>(p.count);

  const auto L1 = exp_marginals1[i1];
  const auto L2 = exp_marginals2[i2];

  const auto exp = std::max(_exp->at(i1, i2), cutoff);

  const auto delta = intra_matrix ? p.coords.bin2.start() - p.coords.bin1.start() : _min_delta;
  if (delta < _min_delta || delta >= _max_delta) {
    _value = {p, exp, 1.0, 0.0, 0.0};
    return _value;
  }

  const auto odds_ratio = compute_odds_ratio(obs, static_cast<double>(obs_sum), N1, N2);
  const auto omega = intra_matrix ? compute_odds_ratio(exp, exp_sum, L1, L2) : 1;

  if ((L1 - exp) * (L2 - exp) <= cutoff) {
    _value = {p, exp, 1.0, odds_ratio, omega};
    return _value;
  }

  if (!std::isfinite(omega) || omega > odds_ratio) {
    _value = {p, exp, 1.0, odds_ratio, omega};
    return _value;
  }

  const auto pvalue =
      compute_pvalue_nchg(*_buffer, static_cast<std::uint64_t>(obs), static_cast<std::uint64_t>(N1),
                          static_cast<std::uint64_t>(N2), obs_sum, omega);
  _value = {p, exp, pvalue, odds_ratio, omega};

  return _value;
}

template <typename File>
inline auto NCHG<File>::iterator::operator->() const -> const_pointer {
  _value = **this;
  return &_value;
}

template <typename File>
inline auto NCHG<File>::iterator::operator++() -> iterator & {
  while (++_pixel_it != _sentinel_it) {
    const auto bin1_id = _pixel_it->coords.bin1.rel_id();
    const auto bin2_id = _pixel_it->coords.bin2.rel_id();

    const auto bin1_masked = (*_bin_mask1)[bin1_id];
    const auto bin2_masked = (*_bin_mask2)[bin2_id];

    if (!bin1_masked && !bin2_masked) {
      break;
    }
  }
  return *this;
}

template <typename File>
inline auto NCHG<File>::iterator::operator++(int) -> iterator {
  auto it = *this;
  it._buffer = std::make_shared<std::vector<double>>();
  std::ignore = ++(*this);
  return it;
}

template <typename File>
inline auto NCHG<File>::compute(const hictk::GenomicInterval &range, double bad_bin_fraction) const
    -> Stats {
  return compute(range, range, bad_bin_fraction);
}

template <typename File>
inline auto NCHG<File>::compute(const hictk::GenomicInterval &range1,
                                const hictk::GenomicInterval &range2, double bad_bin_fraction) const
    -> Stats {
  if (bad_bin_fraction < 0.0 || bad_bin_fraction > 1.0) {
    throw std::logic_error("bad_bin_fraction should be between 0 and 1");
  }

  const auto &chrom1 = range1.chrom();
  const auto &chrom2 = range2.chrom();

  const auto &obs_marginals1 = _obs_matrix->marginals1();
  const auto &obs_marginals2 = _obs_matrix->marginals2();
  const auto &exp_marginals1 = _exp_matrix->marginals1();
  const auto &exp_marginals2 = _exp_matrix->marginals2();

  const auto intra_matrix = chrom1 == chrom2;

  const auto obs_sum = _obs_matrix->sum();
  const auto exp_sum = _exp_matrix->sum();

  const auto &mask1 = *_expected_values.bin_mask(chrom1, chrom2).first;
  const auto &mask2 = *_expected_values.bin_mask(chrom1, chrom2).second;

  const double cutoff = 1.0e-20;

  double N1 = 0.0;
  double N2 = 0.0;
  double L1 = 0.0;
  double L2 = 0.0;

  const auto resolution = _fp->resolution();

  const auto i11 = range1.start() / resolution;
  const auto i12 = (range1.end() + resolution - 1) / resolution;
  const auto i21 = range2.start() / resolution;
  const auto i22 = (range2.end() + resolution - 1) / resolution;

  std::size_t bin1_masked = 0;
  for (auto i = i11; i < i12; ++i) {
    N1 += static_cast<double>(obs_marginals1[i]);
    L1 += exp_marginals1[i];
    bin1_masked += mask1[i];
  }

  std::size_t bin2_masked = 0;
  for (auto i = i21; i < i22; ++i) {
    N2 += static_cast<double>(obs_marginals2[i]);
    L2 += exp_marginals2[i];
    bin2_masked += mask2[i];
  }

  const auto bin1_masked_frac = static_cast<double>(bin1_masked) / static_cast<double>(i12 - i11);
  const auto bin2_masked_frac = static_cast<double>(bin2_masked) / static_cast<double>(i22 - i21);

  double obs = 0.0;
  double exp = 0.0;

  if (bin1_masked_frac < bad_bin_fraction && bin2_masked_frac < bad_bin_fraction) {
    const auto sel = _fp->fetch(range1.chrom().name(), range1.start(), range1.end(),
                                range2.chrom().name(), range2.start(), range2.end());
    const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<double>(),
                                                      sel.template end<double>(), _fp->bins_ptr());

    std::for_each(jsel.begin(), jsel.end(), [&](const hictk::Pixel<double> &p) {
      const auto delta = intra_matrix ? p.coords.bin2.start() - p.coords.bin1.start() : min_delta();

      const auto bin1_id = p.coords.bin1.rel_id();
      const auto bin2_id = p.coords.bin2.rel_id();
      if (delta >= min_delta() && delta < max_delta() && !mask1[bin1_id] && !mask2[bin2_id]) {
        obs += p.count;
        exp += _exp_matrix->at(bin1_id, bin2_id);
      }
    });
  } else {
    obs = 0;
    exp = 0;
  }

  // clang-format off
  const hictk::Pixel<std::uint32_t> p{
      range1.chrom(), range1.start(), range1.end(),
      range2.chrom(), range2.start(), range2.end(),
      static_cast<std::uint32_t>(obs)};
  // clang-format on

  const auto odds_ratio = compute_odds_ratio(obs, static_cast<double>(obs_sum), N1, N2);
  const auto omega = intra_matrix ? compute_odds_ratio(exp, exp_sum, L1, L2) : 1;

  if ((L1 - exp) * (L2 - exp) <= cutoff) {
    return {p, exp, 1.0, odds_ratio, omega};
  }

  if (!std::isfinite(omega) || omega > odds_ratio) {
    return {p, exp, 1.0, odds_ratio, omega};
  }

  const auto pvalue = compute_pvalue_nchg(_nchg_pval_buffer, static_cast<std::uint64_t>(obs),
                                          static_cast<std::uint64_t>(N1),
                                          static_cast<std::uint64_t>(N2), obs_sum, omega);
  return {p, exp, pvalue, odds_ratio, omega};
}

template <typename File>
inline auto NCHG<File>::cbegin(const hictk::Chromosome &chrom1,
                               const hictk::Chromosome &chrom2) const -> iterator {
  const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
  const hictk::transformers::JoinGenomicCoords jsel(
      sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), _fp->bins_ptr());

  return {jsel.begin(),
          jsel.end(),
          _obs_matrix,
          _exp_matrix,
          _expected_values.bin_mask(chrom1),
          _expected_values.bin_mask(chrom2),
          min_delta(),
          max_delta()};
}

template <typename File>
inline auto NCHG<File>::cend(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> iterator {
  const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
  const hictk::transformers::JoinGenomicCoords jsel(
      sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), _fp->bins_ptr());

  return {jsel.end(), jsel.end(), _obs_matrix, _exp_matrix,
          nullptr,    nullptr,    min_delta(), max_delta()};
}

template <typename File>
inline auto NCHG<File>::begin(const hictk::Chromosome &chrom1,
                              const hictk::Chromosome &chrom2) const -> iterator {
  return cbegin(chrom1, chrom2);
}

template <typename File>
inline auto NCHG<File>::end(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> iterator {
  return cend(chrom1, chrom2);
}

template <typename File>
inline auto NCHG<File>::init_exp_matrix(const hictk::Chromosome &chrom1,
                                        const hictk::Chromosome &chrom2, const File &fp,
                                        const ExpectedValues<File> &expected_values)
    -> std::shared_ptr<const ExpectedMatrix<PixelIt>> {
  SPDLOG_INFO(FMT_STRING("[{}:{}] initializing expected matrix..."), chrom1.name(), chrom2.name());

  const auto sel = fp.fetch(chrom1.name(), chrom2.name());
  const auto jsel = hictk::transformers::JoinGenomicCoords(sel.template begin<N>(),
                                                           sel.template end<N>(), fp.bins_ptr());
  const auto first_pixel = jsel.begin();
  const auto last_pixel = jsel.end();

  return std::make_shared<const ExpectedMatrix<PixelIt>>(
      expected_values.expected_matrix(chrom1, chrom2, fp.bins(), first_pixel, last_pixel));
}

template <typename File>
inline auto NCHG<File>::init_obs_matrix(const hictk::Chromosome &chrom1,
                                        const hictk::Chromosome &chrom2, const File &fp,
                                        const std::vector<bool> &bin1_mask,
                                        const std::vector<bool> &bin2_mask, double mad_max_,
                                        std::uint64_t min_delta_, std::uint64_t max_delta_)
    -> std::shared_ptr<const ObservedMatrix<PixelIt>> {
  SPDLOG_INFO(FMT_STRING("[{}:{}] initializing observed matrix..."), chrom1.name(), chrom2.name());

  const auto sel = fp.fetch(chrom1.name(), chrom2.name());
  const auto jsel = hictk::transformers::JoinGenomicCoords(sel.template begin<N>(),
                                                           sel.template end<N>(), fp.bins_ptr());
  const auto first_pixel = jsel.begin();
  const auto last_pixel = jsel.end();

  return std::make_shared<const ObservedMatrix<PixelIt>>(first_pixel, last_pixel, chrom1, chrom2,
                                                         fp.bins(), mad_max_, bin1_mask, bin2_mask,
                                                         min_delta_, max_delta_);
}

template <typename File>
inline double NCHG<File>::compute_cumulative_nchg(std::vector<double> &buffer, std::uint64_t obs,
                                                  std::uint64_t N1, std::uint64_t N2,
                                                  std::uint64_t N, double odds, double precision,
                                                  bool lower_tail) {
  assert(N1 >= obs);
  assert(N2 >= obs);
  assert(N >= obs);
  assert(odds >= 0);
  assert(precision >= 0);

  struct CFishersNCHypergeometricMakeTableParams {
    std::int32_t x1{};
    std::int32_t x2{};
  };

  const auto ub = std::numeric_limits<std::int32_t>::max();

  if (N2 > ub || N1 > ub || N > ub) {
    throw std::runtime_error(
        "unable to compute NCHG CDF: one of the parameters is too large to fit in a "
        "std::int32_t");
  }

  CFishersNCHypergeometric nchg{static_cast<std::int32_t>(N2), static_cast<std::int32_t>(N1),
                                static_cast<std::int32_t>(N), odds, precision};

  CFishersNCHypergeometricMakeTableParams params{};

  const auto size =
      static_cast<std::size_t>(nchg.MakeTable(nullptr, 0, nullptr, nullptr, precision * 0.001));
  buffer.resize(size, 0);

  const auto factor = 1.0 / nchg.MakeTable(buffer.data(), static_cast<std::int32_t>(buffer.size()),
                                           &params.x1, &params.x2, precision * 0.001);
  const auto x_mean = static_cast<std::uint64_t>(lround(nchg.mean()));
  const auto x1 = static_cast<std::uint32_t>(params.x1);
  const auto x2 = static_cast<std::uint32_t>(params.x2);

  assert(x_mean >= x1);
  assert(x_mean <= x2);

  double sum{};
  // Make left tail of table cumulative
  for (std::size_t i = x1; i <= x_mean; ++i) {
    sum = buffer[i - x1] += sum;
  }

  sum = 0.0;
  // Make right tail of table cumulative from the right
  for (std::size_t i = x2; i > x_mean; --i) {
    sum = buffer[i - x1] += sum;
  }

  if (obs <= x_mean) {  // Left tail
    // return 0 when value is outside of table
    const auto pval = obs < x1 ? 0.0 : buffer[obs - x1] * factor;
    return lower_tail ? pval : 1.0 - pval;
  }

  // Right tail, return 0 when value is outside of table
  const auto pval = obs >= x2 ? 0.0 : buffer[obs - x1 + 1] * factor;
  return lower_tail ? 1.0 - pval : pval;
}

template <typename File>
inline double NCHG<File>::compute_pvalue_nchg(std::uint64_t obs, std::uint64_t N1, std::uint64_t N2,
                                              std::uint64_t N, double odds, double precision,
                                              double min_omega) {
  std::vector<double> buffer{};
  return compute_pvalue_nchg(buffer, obs, N1, N2, N, odds, precision, min_omega);
}

template <typename File>
inline double NCHG<File>::compute_pvalue_nchg(std::vector<double> &buffer, std::uint64_t obs,
                                              std::uint64_t N1, std::uint64_t N2, std::uint64_t N,
                                              double odds, double precision, double min_omega) {
  if (obs == 1) {
    return 1.0;
  }
  const auto pvalue =
      compute_cumulative_nchg(buffer, obs - 1, N1, N2, N, std::max(odds, min_omega), precision);
  return pvalue < 0 ? 1.0 : pvalue;
}

template <typename File>
inline auto NCHG<File>::compute_expected_profile() const
    -> std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>> {
  SPDLOG_INFO(FMT_STRING("initializing expected matrix weights from genome-wide interactions..."));
  std::vector<ThinPixelIt> heads{};
  std::vector<ThinPixelIt> tails{};

  using PixelSel = decltype(_fp->fetch("chr1"));
  std::vector<PixelSel> sels{};

  for (const auto &chrom : _fp->chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    auto sel = _fp->fetch(chrom.name());

    auto first = sel.template begin<N>();
    auto last = sel.template end<N>();

    if (first != last) {
      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
      sels.emplace_back(std::move(sel));
    }
  }

  if (heads.empty()) {
    std::vector<double> weights(_fp->bins().size(), 0);
    phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
    for (const auto &chrom : _fp->chromosomes()) {
      scaling_factors.emplace(chrom, 1.0);
    }
    return std::make_pair(weights, scaling_factors);
  }

  hictk::transformers::PixelMerger merger(std::move(heads), std::move(tails));
  const hictk::transformers::JoinGenomicCoords mjsel(merger.begin(), merger.end(), _fp->bins_ptr());

  using PixelItMerged = decltype(mjsel.begin());
  return ExpectedMatrix<PixelItMerged>::build_expected_vector(
      mjsel.begin(), mjsel.end(), _fp->bins(), min_delta(), max_delta());
}

}  // namespace nchg
