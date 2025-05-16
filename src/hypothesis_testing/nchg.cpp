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

#include "nchg/nchg.hpp"

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>
#include <stocc.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/pixel.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <limits>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {

bool NCHGResult::operator<(const NCHGResult &other) const noexcept {
  return pixel.coords < other.pixel.coords;
}

bool NCHGResult::operator>(const NCHGResult &other) const noexcept {
  return pixel.coords > other.pixel.coords;
}

bool NCHGResult::operator==(const NCHGResult &other) const noexcept {
  return pixel.coords == other.pixel.coords;
}

bool NCHGResult::operator!=(const NCHGResult &other) const noexcept { return !(*this == other); }

[[nodiscard]] static std::pair<std::uint64_t, double> aggregate_pixels_cis(
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

[[nodiscard]] static std::pair<std::uint64_t, double> aggregate_pixels_trans(
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

auto NCHG::aggregate_pixels(const hictk::GenomicInterval &range1,
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
    const auto [obs, exp] = aggregate_pixels_cis(*_fp, *_exp_matrix, range1, range2,
                                                 params().min_delta, params().max_delta, *mask);

    return Result{obs, exp};
  }

  const auto &[mask1, mask2] = _expected_values.bin_mask(range1.chrom(), range2.chrom());

  assert(mask1);
  assert(mask2);

  const auto [obs, exp] =
      aggregate_pixels_trans(*_fp, *_exp_matrix, range1, range2, *mask1, *mask2);

  return Result{obs, exp};
}

NCHG::NCHG(std::shared_ptr<const hictk::File> f, const hictk::Chromosome &chrom1,
           const hictk::Chromosome &chrom2, const Params &params)
    : NCHG(f, chrom1, chrom2, ExpectedValues::chromosome_pair(f, chrom1, chrom2, params)) {}

NCHG::NCHG(std::shared_ptr<const hictk::File> f, const hictk::Chromosome &chrom1,
           const hictk::Chromosome &chrom2, ExpectedValues expected_values)
    : _fp(std::move(f)),
      _chrom1(chrom1),
      _chrom2(chrom2),
      _expected_values(std::move(expected_values)) {
  const auto &[bin1_mask, bin2_mask] = _expected_values.bin_mask(chrom1, chrom2);
  auto [obs_matrix, exp_matrix] = init_matrices(
      chrom1, chrom2, *_fp, _expected_values, bin1_mask ? *bin1_mask : std::vector<bool>{},
      bin2_mask ? *bin2_mask : std::vector<bool>{}, _expected_values.params().mad_max,
      _expected_values.params().min_delta, _expected_values.params().max_delta);

  _obs_matrix = std::move(obs_matrix);
  _exp_matrix = std::move(exp_matrix);
}

auto NCHG::params() const noexcept -> Params { return _expected_values.params(); }

const ObservedMatrix &NCHG::observed_matrix() const noexcept { return *_obs_matrix; }

const ExpectedMatrixStats &NCHG::expected_matrix() const noexcept { return *_exp_matrix; }

auto NCHG::compute(const BEDPE &domain, std::uint64_t obs, double exp,
                   double bad_bin_fraction) const -> Stats {
  if (bad_bin_fraction < 0.0 || bad_bin_fraction > 1.0) [[unlikely]] {
    throw std::logic_error("bad_bin_fraction should be between 0 and 1");
  }

  const auto &range1 = domain.range1();
  const auto &range2 = domain.range2();

  const auto obs_sum = _obs_matrix->sum();
  const auto exp_sum = _exp_matrix->sum();

  const auto N1 = compute_N1(range1, range2, bad_bin_fraction);
  const auto N2 = compute_N2(range1, range2, bad_bin_fraction);
  const auto L1 = compute_L1(range1, range2, bad_bin_fraction);
  const auto L2 = compute_L2(range1, range2, bad_bin_fraction);

  const auto range1_masked = N1 == std::numeric_limits<decltype(N1)>::max();
  const auto range2_masked = N2 == std::numeric_limits<decltype(N2)>::max();

  if (range1_masked || range2_masked) [[unlikely]] {
    obs = 0;
    exp = 0.0;
  }

  // clang-format off
  hictk::Pixel p{
      range1.chrom(), range1.start(), range1.end(),
      range2.chrom(), range2.start(), range2.end(),
      obs};
  // clang-format on

  return compute_stats(std::move(p), exp, obs_sum, N1, N2, exp_sum, L1, L2, _nchg_pval_buffer);
}

auto NCHG::cbegin(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> IteratorVariant {
  return std::visit(
      [&](const auto &f) -> IteratorVariant {
        auto [bin1_mask, bin2_mask] = _expected_values.bin_mask(chrom1, chrom2);
        return {iterator{f.fetch(chrom1.name(), chrom2.name()), _obs_matrix, _exp_matrix, bin1_mask,
                         bin2_mask, params().min_delta, params().max_delta}};
      },
      _fp->get());
}

auto NCHG::cend(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> IteratorVariant {
  return std::visit(
      [&](const auto &f) -> IteratorVariant {
        using PixelSelector = decltype(f.fetch("chr1", "chr2"));
        return {iterator<PixelSelector>::at_end(f.fetch(chrom1.name(), chrom2.name()), _obs_matrix,
                                                _exp_matrix)};
      },
      _fp->get());
}

auto NCHG::begin(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> IteratorVariant {
  return cbegin(chrom1, chrom2);
}

auto NCHG::end(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> IteratorVariant {
  return cend(chrom1, chrom2);
}

auto NCHG::init_matrices(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
                         const hictk::File &f, const ExpectedValues &expected_values,
                         const std::vector<bool> &bin1_mask, const std::vector<bool> &bin2_mask,
                         double mad_max_, std::uint64_t min_delta_, std::uint64_t max_delta_)
    -> std::pair<std::shared_ptr<const ObservedMatrix>,
                 std::shared_ptr<const ExpectedMatrixStats>> {
  SPDLOG_INFO("[{}:{}] initializing observed and expected matrices...", chrom1.name(),
              chrom2.name());

  MatrixStats<std::uint32_t> obs_stats{_chrom1,        _chrom2,    bin1_mask, bin2_mask,
                                       f.resolution(), min_delta_, max_delta_};
  MatrixStats<double> exp_stats{_chrom1,        _chrom2,    bin1_mask,  bin2_mask,
                                f.resolution(), min_delta_, max_delta_, expected_values.weights()};

  std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom1.name(), chrom2.name());
        const auto jsel = hictk::transformers::JoinGenomicCoords(
            sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), fp.bins_ptr());

        for (const auto &p : jsel) {
          obs_stats.add(p);
          exp_stats.add(hictk::Pixel{p.coords, static_cast<double>(p.count)});
        }
      },
      f.get());

  auto obs_matrix = std::make_shared<const ObservedMatrix>(
      chrom1, chrom2, f.bins(), std::move(obs_stats.marginals1), std::move(obs_stats.marginals2),
      obs_stats.nnz, obs_stats.sum, mad_max_, min_delta_, max_delta_);

  auto exp_matrix = std::make_shared<const ExpectedMatrixStats>(
      chrom1, chrom2, f.bins(), expected_values.weights(), std::move(exp_stats.marginals1),
      std::move(exp_stats.marginals2), exp_stats.nnz, exp_stats.sum, min_delta_, max_delta_);

  SPDLOG_INFO(
      "[{}:{}] initialed observed and expected matrices with a total of {} interactions ({} nnz)!",
      chrom1.name(), chrom2.name(), obs_stats.sum, obs_stats.nnz);

  return {std::move(obs_matrix), std::move(exp_matrix)};
}

double NCHG::compute_cumulative_nchg(std::vector<double> &buffer, std::uint64_t obs,
                                     std::uint64_t N1, std::uint64_t N2, std::uint64_t N,
                                     double odds, double precision, bool lower_tail) {
  assert(N1 >= obs);
  assert(N2 >= obs);
  assert(N >= obs);
  assert(odds >= 0);
  assert(precision >= 0);

  struct CFishersNCHypergeometricMakeTableParams {
    std::int32_t x1{};
    std::int32_t x2{};
  };

  constexpr auto ub = std::numeric_limits<std::int32_t>::max();

  if (N2 > ub || N1 > ub || N > ub) [[unlikely]] {
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
  const auto x_mean = static_cast<std::uint64_t>(std::lround(nchg.mean()));
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

double NCHG::compute_pvalue_nchg(std::uint64_t obs, std::uint64_t N1, std::uint64_t N2,
                                 std::uint64_t N, double odds, double precision, double min_omega) {
  std::vector<double> buffer{};
  return compute_pvalue_nchg(buffer, obs, N1, N2, N, odds, precision, min_omega);
}

double NCHG::compute_pvalue_nchg(std::vector<double> &buffer, std::uint64_t obs, std::uint64_t N1,
                                 std::uint64_t N2, std::uint64_t N, double odds, double precision,
                                 double min_omega) {
  if (obs == 1) [[unlikely]] {
    return 1.0;
  }
  const auto pvalue =
      compute_cumulative_nchg(buffer, obs - 1, N1, N2, N, std::max(odds, min_omega), precision);
  return pvalue < 0 ? 1.0 : pvalue;
}

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
    return {0.0, N2{}};
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

std::uint64_t NCHG::compute_N1(const hictk::GenomicInterval &range1,
                               const hictk::GenomicInterval &range2,
                               double max_bad_bin_threshold) const noexcept {
  assert(_obs_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] =
      aggregate_marginals(range1, _fp->resolution(), _obs_matrix->marginals1(),
                          *_expected_values.bin_mask(range1.chrom(), range2.chrom()).first);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<std::uint64_t>::max();
  }

  return sum;
}

std::uint64_t NCHG::compute_N2(const hictk::GenomicInterval &range1,
                               const hictk::GenomicInterval &range2,
                               double max_bad_bin_threshold) const noexcept {
  assert(_obs_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] =
      aggregate_marginals(range2, _fp->resolution(), _obs_matrix->marginals2(),
                          *_expected_values.bin_mask(range1.chrom(), range2.chrom()).second);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<std::uint64_t>::max();
  }

  return sum;
}

double NCHG::compute_L1(const hictk::GenomicInterval &range1, const hictk::GenomicInterval &range2,
                        double max_bad_bin_threshold) const noexcept {
  assert(_exp_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] =
      aggregate_marginals(range1, _fp->resolution(), _exp_matrix->marginals1(),
                          *_expected_values.bin_mask(range1.chrom(), range2.chrom()).first);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return sum;
}

double NCHG::compute_L2(const hictk::GenomicInterval &range1, const hictk::GenomicInterval &range2,
                        double max_bad_bin_threshold) const noexcept {
  assert(_exp_matrix);
  assert(_fp);
  const auto &[bad_bin_frac, sum] =
      aggregate_marginals(range2, _fp->resolution(), _exp_matrix->marginals2(),
                          *_expected_values.bin_mask(range1.chrom(), range2.chrom()).second);

  if (bad_bin_frac >= max_bad_bin_threshold) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return sum;
}

NCHGResult NCHG::compute_stats(hictk::Pixel<std::uint64_t> pixel, double exp, std::uint64_t obs_sum,
                               std::uint64_t N1, std::uint64_t N2, double exp_sum, double L1,
                               double L2, std::vector<double> &buffer) {
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

  const auto odds_ratio = compute_odds_ratio(pixel.count, obs_sum, N1, N2);
  const auto omega = intra_matrix ? compute_odds_ratio(exp, exp_sum, L1, L2) : 1.0;
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

  const auto pvalue = compute_pvalue_nchg(buffer, pixel.count, N1, N2, obs_sum, omega);
  return {.pixel = std::move(pixel),
          .expected = exp,
          .pval = pvalue,
          .log_ratio = log_ratio,
          .odds_ratio = odds_ratio,
          .omega = omega};
}

}  // namespace nchg
