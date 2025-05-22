// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>

#pragma once

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace nchg {

class ExpectedValuesAggregator {
  std::shared_ptr<const hictk::BinTable> _bins;
  std::size_t _num_bins_gw{};

  using CisKey = hictk::Chromosome;
  using TransKey = std::pair<CisKey, CisKey>;
  phmap::flat_hash_map<CisKey, double> _cis_sum;
  phmap::flat_hash_map<TransKey, double> _trans_sum;

  std::vector<double> _possible_distances;
  std::vector<double> _observed_distances;

  std::vector<double> _weights;
  phmap::btree_map<hictk::Chromosome, double> _scaling_factors;

 public:
  ExpectedValuesAggregator() = default;
  explicit ExpectedValuesAggregator(std::shared_ptr<const hictk::BinTable> bins);
  template <typename N>
  void add(const hictk::ThinPixel<N>& p);
  template <typename N>
  void add(const hictk::Pixel<N>& p);

  void compute_density(double bin_aggregation_possible_distances_cutoff,
                       double bin_aggregation_observed_distances_cutoff, bool intepolate_weights,
                       double interpolation_qtile, std::uint32_t interpolation_window_size);

  [[nodiscard]] const std::vector<double>& weights() const noexcept;
  [[nodiscard]] std::vector<double> weights(const hictk::Chromosome& chrom,
                                            bool rescale = true) const;

  [[nodiscard]] double scaling_factor(const hictk::Chromosome& chrom) const;
  [[nodiscard]] const phmap::btree_map<hictk::Chromosome, double>& scaling_factors() const noexcept;

 private:
  [[nodiscard]] const hictk::Reference& chromosomes() const noexcept;

  void init_possible_distances();
  void compute_density_cis(double bin_aggregation_possible_distances_cutoff,
                           double bin_aggregation_observed_distances_cutoff,
                           bool interpolate_weights, double interpolation_qtile,
                           std::uint32_t interpolation_window_size);
  void compute_density_trans();

  std::size_t aggregate_low_cov_bins(double possible_distances_cutoff,
                                     double observed_distances_cutoff) noexcept;
  [[nodiscard]] phmap::btree_map<hictk::Chromosome, double> compute_scaling_factors() const;
  void correct_outliers(double quantile, std::uint32_t window_size);

  [[nodiscard]] double at(const hictk::Chromosome& chrom) const;
  [[nodiscard]] double at(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2) const;

  [[nodiscard]] double& at(const hictk::Chromosome& chrom);
  [[nodiscard]] double& at(const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2);
};

}  // namespace nchg

#include "../../expected_values_aggregator_impl.hpp"  // NOLINT
