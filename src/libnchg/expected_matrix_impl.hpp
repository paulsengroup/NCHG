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

#pragma once

#include <algorithm>
#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/hic/expected_values_aggregator.hpp>
#include <hictk/pixel.hpp>
#include <utility>
#include <variant>
#include <vector>

namespace nchg {

template <typename PixelIt>
inline auto ExpectedMatrix<PixelIt>::compute_stats(PixelIt first_pixel, PixelIt last_pixel,
                                                   const hictk::Chromosome &chrom1,
                                                   const hictk::Chromosome &chrom2,
                                                   const hictk::BinTable &bins,
                                                   const std::vector<double> &weights,
                                                   double scaling_factor, std::uint64_t min_delta_,
                                                   std::uint64_t max_delta_) {
  struct Result {
    std::shared_ptr<std::vector<double>> marginals1{std::make_shared<std::vector<double>>()};
    std::shared_ptr<std::vector<double>> marginals2{std::make_shared<std::vector<double>>()};
    double sum{};
    std::uint64_t nnz{};
  };

  Result res{};
  const auto intra_matrix = chrom1 == chrom2;
  if (intra_matrix) {
    res.marginals2 = res.marginals1;
  }

  const auto bin_size = bins.resolution();
  const auto num_bins1 = (chrom1.size() + bin_size - 1) / bin_size;
  const auto num_bins2 = (chrom2.size() + bin_size - 1) / bin_size;

  if (num_bins1 == 0 || num_bins2 == 0) {
    return res;
  }

  auto &marginals1_ = *res.marginals1;
  auto &marginals2_ = *res.marginals2;
  auto &sum_ = res.sum;
  auto &nnz_ = res.nnz;
  marginals1_.resize(num_bins1, 0);
  marginals2_.resize(num_bins2, 0);

  std::for_each(first_pixel, last_pixel, [&](const hictk::Pixel<N> &p) {
    const auto &bin1 = p.coords.bin1;
    const auto &bin2 = p.coords.bin2;
    const auto delta = bin2.start() - bin1.start();
    if (intra_matrix && (delta < min_delta_ || delta >= max_delta_)) {
      return;
    }

    auto count = intra_matrix ? weights.at(bin2.id() - bin1.id()) / scaling_factor
                              : static_cast<double>(p.count);
    if (std::isnan(count)) {
      count = 0.0;
    }

    if (intra_matrix) {
      sum_ += bin1 == bin2 ? count : 2 * count;
      nnz_ += bin1 == bin2 ? 1ULL : 2ULL;
    } else {
      sum_ += count;
      nnz_++;
    }

    const auto i1 = bin1.rel_id();
    marginals1_[i1] += count;

    if (!intra_matrix || bin1 != bin2) {
      const auto i2 = bin2.rel_id();
      marginals2_[i2] += count;
    }
  });

  return res;
}

template <typename PixelIt>
template <typename PixelItGw>
inline ExpectedMatrix<PixelIt>::ExpectedMatrix(PixelIt first_pixel, PixelIt last_pixel,
                                               PixelItGw first_pixel_gw, PixelItGw last_pixel_gw,
                                               const hictk::Chromosome &chrom1,
                                               const hictk::Chromosome &chrom2,
                                               const hictk::BinTable &bins,
                                               std::uint64_t min_delta_, std::uint64_t max_delta_)
    : ExpectedMatrix(std::move(first_pixel), std::move(last_pixel), chrom1, chrom2, bins,
                     compute_weights(std::move(first_pixel_gw), std::move(last_pixel_gw), chrom1,
                                     chrom2, bins, min_delta_, max_delta_),
                     min_delta_, max_delta_) {}

template <typename PixelIt>
inline ExpectedMatrix<PixelIt>::ExpectedMatrix(PixelIt first_pixel, PixelIt last_pixel,
                                               hictk::Chromosome chrom1, hictk::Chromosome chrom2,
                                               hictk::BinTable bins, std::vector<double> weights,
                                               double scaling_factor, std::uint64_t min_delta_,
                                               std::uint64_t max_delta_)
    : _first(std::move(first_pixel)),
      _last(std::move(last_pixel)),
      _chrom1(std::move(chrom1)),
      _chrom2(std::move(chrom2)),
      _bins(std::move(bins)),
      _weights(std::move(weights)),
      _min_delta(min_delta_),
      _max_delta(max_delta_) {
  auto stats = compute_stats(_first, _last, _chrom1, _chrom2, _bins, _weights, scaling_factor,
                             min_delta_, max_delta_);
  _marginals1 = std::move(stats.marginals1);
  _marginals2 = std::move(stats.marginals2);
  _nnz = stats.nnz;
  _sum = stats.sum;
}

template <typename PixelIt>
inline std::uint32_t ExpectedMatrix<PixelIt>::resolution() const noexcept {
  return _bins.resolution();
}

template <typename PixelIt>
inline std::size_t ExpectedMatrix<PixelIt>::num_rows() const noexcept {
  return _bins.subset(_chrom1).size();
}

template <typename PixelIt>
inline std::size_t ExpectedMatrix<PixelIt>::num_cols() const noexcept {
  return _bins.subset(_chrom2).size();
}

template <typename PixelIt>
inline const hictk::Chromosome &ExpectedMatrix<PixelIt>::chrom1() const noexcept {
  return _chrom1;
}
template <typename PixelIt>
inline const hictk::Chromosome &ExpectedMatrix<PixelIt>::chrom2() const noexcept {
  return _chrom2;
}

template <typename PixelIt>
inline std::uint64_t ExpectedMatrix<PixelIt>::nnz() const noexcept {
  return _nnz;
}

template <typename PixelIt>
inline double ExpectedMatrix<PixelIt>::sum() const noexcept {
  return _sum;
}

template <typename PixelIt>
inline double ExpectedMatrix<PixelIt>::nnz_avg() const noexcept {
  return sum() / static_cast<double>(nnz());
}

template <typename PixelIt>
inline const std::vector<double> &ExpectedMatrix<PixelIt>::weights() const noexcept {
  return _weights;
}

template <typename PixelIt>
inline const phmap::btree_map<hictk::Chromosome, double> &ExpectedMatrix<PixelIt>::scaling_factors()
    const noexcept {
  return _scaling_factors;
}

template <typename PixelIt>
inline std::uint64_t ExpectedMatrix<PixelIt>::min_delta() const noexcept {
  return _min_delta;
}

template <typename PixelIt>
inline std::uint64_t ExpectedMatrix<PixelIt>::max_delta() const noexcept {
  return _max_delta;
}

template <typename PixelIt>
inline double ExpectedMatrix<PixelIt>::at(std::uint64_t i, std::uint64_t j) const {
  if (_chrom1 == _chrom2) {
    return _weights.at(j - i);
  }
  return nnz_avg();
}

template <typename PixelIt>
inline const std::vector<double> &ExpectedMatrix<PixelIt>::marginals1() const noexcept {
  return *_marginals1;
}
template <typename PixelIt>
inline const std::vector<double> &ExpectedMatrix<PixelIt>::marginals2() const noexcept {
  return *_marginals2;
}

template <typename PixelIt>
template <typename PixelItGw>
inline std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
ExpectedMatrix<PixelIt>::build_expected_vector(PixelItGw first_pixel, PixelItGw last_pixel,
                                               const hictk::Chromosome &chrom1,
                                               const hictk::Chromosome &chrom2,
                                               const hictk::BinTable &bins,
                                               std::uint64_t min_delta_, std::uint64_t max_delta_) {
  if (chrom1 != chrom2) {
    return {};
  }

  SPDLOG_INFO(FMT_STRING("building expected value vector using hictk..."));

  if (first_pixel == last_pixel) {
    phmap::btree_map<hictk::Chromosome, double> scaling_factors{};
    std::transform(bins.chromosomes().begin(), bins.chromosomes().end(),
                   std::inserter(scaling_factors, scaling_factors.begin()),
                   [](const hictk::Chromosome &chrom) { return std::make_pair(chrom, 0.0); });
    return std::make_pair(std::vector<double>(bins.size(), 0), scaling_factors);
  }

  const auto bins_ = std::make_shared<const hictk::BinTable>(bins);

  hictk::hic::internal::ExpectedValuesAggregator aggr(bins_);

  std::for_each(first_pixel, last_pixel, [&](const hictk::Pixel<N> &p) {
    const auto delta = p.coords.bin2.start() - p.coords.bin1.start();
    if (delta >= min_delta_ && delta < max_delta_) {
      // TODO remove redundant pixel copy
      hictk::Pixel<float> pp(*bins_, hictk::ThinPixel<float>{p.coords.bin1.id(), p.coords.bin2.id(),
                                                             static_cast<float>(p.count)});
      aggr.add(std::move(pp));
    }
  });
  aggr.compute_density();

  auto weights = aggr.weights();
  const auto &chrom = bins.chromosomes().longest_chromosome();
  weights.resize((chrom.size() + bins.resolution() - 1) / bins.resolution(),
                 std::numeric_limits<double>::quiet_NaN());

  return std::make_pair(std::move(weights), aggr.scaling_factors());
}

template <typename PixelIt>
template <typename PixelItGw>
inline std::vector<double> ExpectedMatrix<PixelIt>::compute_weights(
    PixelItGw first_pixel, PixelItGw last_pixel, const hictk::Chromosome &chrom1,
    const hictk::Chromosome &chrom2, const hictk::BinTable &bins, std::uint64_t min_delta_,
    std::uint64_t max_delta_) {
  if (chrom1 != chrom2) {
    return {};
  }

  auto [weights, scaling_factors] = build_expected_vector(
      std::move(first_pixel), std::move(last_pixel), chrom1, chrom2, bins, min_delta_, max_delta_);

  weights.resize((chrom1.size() + bins.resolution() - 1) / bins.resolution(),
                 std::numeric_limits<double>::quiet_NaN());

  const auto sf = scaling_factors.at(chrom1);
  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [&](const auto n) { return n / sf; });

  return weights;
}

}  // namespace nchg
