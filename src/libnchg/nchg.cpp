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
#include <hictk/transformers/join_genomic_coords.hpp>
#include <memory>
#include <ranges>
#include <stdexcept>
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

NCHG::NCHG(std::shared_ptr<const hictk::File> f, const hictk::Chromosome &chrom1,
           const hictk::Chromosome &chrom2, const Params &params)
    : NCHG(f, chrom1, chrom2, ExpectedValues::chromosome_pair(f, chrom1, chrom2, params)) {}

NCHG::NCHG(std::shared_ptr<const hictk::File> f, const hictk::Chromosome &chrom1,
           const hictk::Chromosome &chrom2, ExpectedValues expected_values)
    : _fp(std::move(f)),
      _chrom1(chrom1),
      _chrom2(chrom2),
      _exp_matrix(init_exp_matrix(chrom1, chrom2, *_fp, expected_values)),
      _obs_matrix(init_obs_matrix(
          chrom1, chrom2, *_fp, *expected_values.bin_mask(chrom1, chrom2).first,
          *expected_values.bin_mask(chrom1, chrom2).second, expected_values.params().mad_max,
          expected_values.params().min_delta, expected_values.params().max_delta)),
      _expected_values(std::move(expected_values)) {}

auto NCHG::params() const noexcept -> Params { return _expected_values.params(); }

const ObservedMatrix &NCHG::observed_matrix() const noexcept { return *_obs_matrix; }

const ExpectedMatrixStats &NCHG::expected_matrix() const noexcept { return *_exp_matrix; }

auto NCHG::compute(const hictk::GenomicInterval &range, double bad_bin_fraction) const -> Stats {
  return compute(range, range, bad_bin_fraction);
}

auto NCHG::compute(const hictk::GenomicInterval &range1, const hictk::GenomicInterval &range2,
                   double bad_bin_fraction) const -> Stats {
  if (bad_bin_fraction < 0.0 || bad_bin_fraction > 1.0) [[unlikely]] {
    throw std::logic_error("bad_bin_fraction should be between 0 and 1");
  }

  const auto obs_sum = _obs_matrix->sum();
  const auto exp_sum = _exp_matrix->sum();

  const auto N1 = compute_N1(range1, range2, bad_bin_fraction);
  const auto N2 = compute_N2(range1, range2, bad_bin_fraction);
  const auto L1 = compute_L1(range1, range2, bad_bin_fraction);
  const auto L2 = compute_L2(range1, range2, bad_bin_fraction);

  const auto range1_masked = std::isnan(N1);
  const auto range2_masked = std::isnan(N1);

  std::uint64_t obs = 0;
  double exp = 0.0;

  if (!range1_masked && !range2_masked) [[likely]] {
    const auto &result = aggregate_pixels(range1, range2);
    obs = result.obs;
    exp = std::max(0.0, result.exp);
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
        return {iterator{f.fetch(chrom1.name(), chrom2.name()), _obs_matrix, _exp_matrix,
                         _expected_values.bin_mask(chrom1, chrom2).first,
                         _expected_values.bin_mask(chrom1, chrom2).second, params().min_delta,
                         params().max_delta}};
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

auto NCHG::init_exp_matrix(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
                           const hictk::File &f, const ExpectedValues &expected_values)
    -> std::shared_ptr<const ExpectedMatrixStats> {
  SPDLOG_INFO("[{}:{}] initializing expected matrix...", chrom1.name(), chrom2.name());

  return std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom1.name(), chrom2.name());
        const auto jsel = hictk::transformers::JoinGenomicCoords(
            sel.template begin<double>(), sel.template end<double>(), fp.bins_ptr());
        return std::make_shared<const ExpectedMatrixStats>(
            expected_values.expected_matrix(chrom1, chrom2, fp.bins(), jsel));
      },
      f.get());
}

auto NCHG::init_obs_matrix(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2,
                           const hictk::File &f, const std::vector<bool> &bin1_mask,
                           const std::vector<bool> &bin2_mask, double mad_max_,
                           std::uint64_t min_delta_, std::uint64_t max_delta_)
    -> std::shared_ptr<const ObservedMatrix> {
  SPDLOG_INFO("[{}:{}] initializing observed matrix...", chrom1.name(), chrom2.name());

  return std::visit(
      [&](const auto &fp) {
        const auto sel = fp.fetch(chrom1.name(), chrom2.name());
        const auto jsel = hictk::transformers::JoinGenomicCoords(
            sel.template begin<N>(), sel.template end<N>(), fp.bins_ptr());
        return std::make_shared<const ObservedMatrix>(jsel, chrom1, chrom2, fp.bins(), mad_max_,
                                                      bin1_mask, bin2_mask, min_delta_, max_delta_);
      },
      f.get());
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

std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
NCHG::compute_expected_profile() const {
  SPDLOG_INFO("initializing expected matrix weights from genome-wide interactions...");

  return std::visit(
      [&](const auto &f)
          -> std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>> {
        const auto [selectors, merger] = ExpectedValues::init_pixel_merger_cis(f);
        if (merger.begin() == merger.end()) [[unlikely]] {
          phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
          std::ranges::transform(
              f.chromosomes(), std::inserter(scaling_factors, scaling_factors.end()),
              [](const hictk::Chromosome &chrom) { return std::make_pair(chrom, 1.0); });
          for (const auto &chrom : f.chromosomes()) {
            scaling_factors.emplace(chrom, 1.0);
          }
          return {std::vector<double>(f.bins().size(), 0), scaling_factors};
        }

        const hictk::transformers::JoinGenomicCoords mjsel(merger.begin(), merger.end(),
                                                           f.bins_ptr());

        return ExpectedMatrixStats::build_expected_vector(mjsel, f.bins(), params().min_delta,
                                                          params().max_delta);
      },
      _fp->get());
}

}  // namespace nchg
