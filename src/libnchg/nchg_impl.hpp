// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <stocc.h>

#include <cassert>
#include <hictk/file.hpp>
#include <hictk/fmt/genomic_interval.hpp>
#include <hictk/hic/expected_values_aggregator.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <hictk/transformers/pixel_merger.hpp>

#include "nchg/expected_matrix.hpp"
#include "nchg/observed_matrix.hpp"

namespace nchg {

template <typename File>
inline NCHG<File>::NCHG(std::shared_ptr<const File> f, std::uint64_t min_delta,
                        std::uint64_t max_delta)
    : _fp(std::move(f)), _min_delta(min_delta), _max_delta(max_delta), _expected_values(_fp) {
  if (_min_delta >= _max_delta) {
    throw std::logic_error("min_delta should be strictly less than max_delta");
  }
}

template <typename File>
inline NCHG<File> NCHG<File>::cis_only(std::shared_ptr<const File> f, std::uint64_t min_delta,
                                       std::uint64_t max_delta) {
  assert(f);
  if (min_delta >= max_delta) {
    throw std::logic_error("min_delta should be strictly less than max_delta");
  }

  NCHG<File> nchg{nullptr};
  nchg._fp = std::move(f);
  nchg._min_delta = min_delta;
  nchg._max_delta = max_delta;
  nchg._expected_values = ExpectedValues<File>::cis_only(nchg._fp, min_delta, max_delta);

  for (const auto &chrom : nchg._fp->chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    nchg.init_matrix(chrom);
  }

  return nchg;
}

template <typename File>
inline NCHG<File> NCHG<File>::trans_only(std::shared_ptr<const File> f) {
  assert(f);

  NCHG<File> nchg{nullptr};
  nchg._fp = std::move(f);
  nchg._expected_values = ExpectedValues<File>::trans_only(nchg._fp);

  const auto num_chroms = nchg._fp->chromosomes().size();
  for (std::uint32_t chrom1_id = 0; chrom1_id < num_chroms; ++chrom1_id) {
    const auto &chrom1 = nchg._fp->chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < num_chroms; ++chrom2_id) {
      const auto &chrom2 = nchg._fp->chromosomes().at(chrom2_id);
      nchg.init_matrix(chrom1, chrom2);
    }
  }
  return nchg;
}

template <typename File>
inline NCHG<File> NCHG<File>::chromosome_pair(std::shared_ptr<const File> f,
                                              const hictk::Chromosome &chrom1,
                                              const hictk::Chromosome &chrom2,
                                              std::uint64_t min_delta, std::uint64_t max_delta) {
  assert(f);
  if (min_delta >= max_delta) {
    throw std::logic_error("min_delta should be strictly less than max_delta");
  }

  NCHG<File> nchg{nullptr};
  nchg._fp = std::move(f);
  nchg._min_delta = min_delta;
  nchg._max_delta = max_delta;
  nchg._expected_values =
      ExpectedValues<File>::chromosome_pair(nchg._fp, chrom1, chrom2, min_delta, max_delta);

  nchg.init_matrix(chrom1, chrom2);

  return nchg;
}

template <typename File>
inline auto NCHG<File>::observed_matrix(const hictk::Chromosome &chrom) const
    -> const ObservedMatrix<PixelIt> & {
  return observed_matrix(chrom, chrom);
}

template <typename File>
inline auto NCHG<File>::observed_matrix(const hictk::Chromosome &chrom1,
                                        const hictk::Chromosome &chrom2) const
    -> const ObservedMatrix<PixelIt> & {
  return *_obs_matrices->at({chrom1, chrom2});
}

template <typename File>
inline auto NCHG<File>::expected_matrix(const hictk::Chromosome &chrom) const
    -> const ExpectedMatrix<PixelIt> & {
  return expected_matrix(chrom, chrom);
}

template <typename File>
inline auto NCHG<File>::expected_matrix(const hictk::Chromosome &chrom1,
                                        const hictk::Chromosome &chrom2) const
    -> const ExpectedMatrix<PixelIt> & {
  return *_exp_matrices->at({chrom1, chrom2});
}

template <typename File>
inline void NCHG<File>::init_matrices() {
  for (std::uint32_t chrom1_id = 0; chrom1_id < _fp->chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = _fp->chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < _fp->chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = _fp->chromosomes().at(chrom2_id);
      init_matrix(chrom1, chrom2);
    }
  }
}

template <typename File>
inline void NCHG<File>::init_matrix(const hictk::Chromosome &chrom) {
  init_matrix(chrom, chrom);
}

template <typename File>
inline void NCHG<File>::init_matrix(const hictk::Chromosome &chrom1,
                                    const hictk::Chromosome &chrom2) {
  const Key k{chrom1, chrom2};
  if (_obs_matrices.find(k) != _obs_matrices.end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("matrix for \"{}:{}\" has already been initialized"), chrom1.name(),
                    chrom2.name()));
  }

  const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
  const auto jsel = hictk::transformers::JoinGenomicCoords(sel.template begin<N>(),
                                                           sel.template end<N>(), _fp->bins_ptr());
  const auto first_pixel = jsel.begin();
  const auto last_pixel = jsel.end();

  SPDLOG_INFO(FMT_STRING("[{}:{}] initializing observed matrix..."), chrom1.name(), chrom2.name());
  _obs_matrices.emplace(
      k, std::make_shared<const ObservedMatrix<PixelIt>>(first_pixel, last_pixel, chrom1, chrom2,
                                                         _fp->bins(), _min_delta, _max_delta));

  assert(_exp_matrices.find(k) == _exp_matrices.end());

  SPDLOG_INFO(FMT_STRING("[{}:{}] initializing expected matrix..."), chrom1.name(), chrom2.name());

  _exp_matrices.emplace(k, std::make_shared<const ExpectedMatrix<PixelIt>>(
                               _expected_values.expected_matrix(chrom1, chrom2)));
}

template <typename File>
inline void NCHG<File>::erase_matrices() noexcept {
  std::vector<Key> chrom_pairs{};
  std::transform(_obs_matrices.begin(), _obs_matrices.end(), std::back_inserter(chrom_pairs),
                 [](const auto &kv) { return kv.first; });

  for (const auto &chrom_pair : chrom_pairs) {
    erase_matrix(chrom_pair.first, chrom_pair.second);
  }
}

template <typename File>
inline void NCHG<File>::erase_matrix(const hictk::Chromosome &chrom) {
  erase_matrix(chrom, chrom);
}

template <typename File>
inline void NCHG<File>::erase_matrix(const hictk::Chromosome &chrom1,
                                     const hictk::Chromosome &chrom2) {
  const Key k{chrom1, chrom2};

  _obs_matrices.erase(k);
  _exp_matrices.erase(k);
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
inline NCHG<File>::iterator::iterator(PixelIt pixel_it,
                                      std::shared_ptr<const ObservedMatrix<PixelIt>> obs,
                                      std::shared_ptr<const ExpectedMatrix<PixelIt>> exp,
                                      std::uint64_t min_delta, std::uint64_t max_delta)
    : _pixel_it(std::move(pixel_it)),
      _obs(std::move(obs)),
      _exp(std::move(exp)),
      _min_delta(min_delta),
      _max_delta(max_delta) {
  assert(min_delta < max_delta);
}

template <typename File>
inline NCHG<File>::iterator::iterator(const iterator &other)
    : _pixel_it(other._pixel_it),
      _obs(other._obs),
      _exp(other._exp),
      _min_delta(other._min_delta),
      _max_delta(other._max_delta),
      _value((other._value)) {}

template <typename File>
inline auto NCHG<File>::iterator::operator=(const iterator &other) -> iterator & {
  if (this == &other) {
    return *this;
  }

  _pixel_it = other._pixel_it;
  _obs = other._obs;
  _exp = other._exp;
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
  ++_pixel_it;
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
inline auto NCHG<File>::compute(const hictk::GenomicInterval &range) const -> Stats {
  return compute(range, range);
}

template <typename File>
inline auto NCHG<File>::compute(const hictk::GenomicInterval &range1,
                                const hictk::GenomicInterval &range2) const -> Stats {
  const auto &chrom1 = range1.chrom();
  const auto &chrom2 = range2.chrom();

  const auto &obs_matrix = *_obs_matrices.at(Key{chrom1, chrom2});
  const auto &exp_matrix = *_exp_matrices.at(Key{chrom1, chrom2});

  const auto &obs_marginals1 = obs_matrix.marginals1();
  const auto &obs_marginals2 = obs_matrix.marginals2();
  const auto &exp_marginals1 = exp_matrix.marginals1();
  const auto &exp_marginals2 = exp_matrix.marginals2();

  const auto intra_matrix = chrom1 == chrom2;

  const auto obs_sum = obs_matrix.sum();
  const auto exp_sum = exp_matrix.sum();

  const double cutoff = 1.0e-20;

  double N1 = 0.0;
  double N2 = 0.0;
  double L1 = 0.0;
  double L2 = 0.0;

  const auto resolution = _fp->resolution();

  // TODO how do we handle partial overlaps?
  const auto i11 = range1.start() / resolution;
  const auto i12 = (range1.end() + resolution - 1) / resolution;
  const auto i21 = range2.start() / resolution;
  const auto i22 = (range2.end() + resolution - 1) / resolution;

  for (auto i = i11; i < i12; ++i) {
    N1 += static_cast<double>(obs_marginals1[i]);
    L1 += exp_marginals1[i];
  }
  for (auto i = i21; i < i22; ++i) {
    N2 += static_cast<double>(obs_marginals2[i]);
    L2 += exp_marginals2[i];
  }

  const auto sel = _fp->fetch(range1.chrom().name(), range1.start(), range1.end(),
                              range2.chrom().name(), range2.start(), range2.end());
  const hictk::transformers::JoinGenomicCoords jsel(sel.template begin<double>(),
                                                    sel.template end<double>(), _fp->bins_ptr());

  double obs = 0.0;
  double exp = 0.0;

  std::for_each(jsel.begin(), jsel.end(), [&](const hictk::Pixel<double> &p) {
    const auto delta = intra_matrix ? p.coords.bin2.start() - p.coords.bin1.start() : _min_delta;
    if (delta >= _min_delta && delta < _max_delta) {
      obs += p.count;
      exp += exp_matrix.at(p.coords.bin1.rel_id(), p.coords.bin2.rel_id());
    }
  });

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
  const auto &obs = _obs_matrices.at(Key{chrom1, chrom2});
  const auto &exp = _exp_matrices.at(Key{chrom1, chrom2});

  const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
  const hictk::transformers::JoinGenomicCoords jsel(
      sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), _fp->bins_ptr());

  return {jsel.begin(), obs, exp, _min_delta, _max_delta};
}

template <typename File>
inline auto NCHG<File>::cend(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
    -> iterator {
  const auto &obs = _obs_matrices.at(Key{chrom1, chrom2});
  const auto &exp = _exp_matrices.at(Key{chrom1, chrom2});

  const auto sel = _fp->fetch(chrom1.name(), chrom2.name());
  const hictk::transformers::JoinGenomicCoords jsel(
      sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), _fp->bins_ptr());

  return {jsel.end(), obs, exp, _min_delta, _max_delta};
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
  return ExpectedMatrix<PixelItMerged>::build_expected_vector(mjsel.begin(), mjsel.end(),
                                                              _fp->bins(), _min_delta, _max_delta);
}

}  // namespace nchg
