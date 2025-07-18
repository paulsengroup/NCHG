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

#include <parallel_hashmap/phmap.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <type_traits>
#include <utility>
#include <vector>

#include "nchg/concepts.hpp"
#include "nchg/expected_matrix.hpp"
#include "nchg/mad_max_filter.hpp"

namespace nchg {

template <typename Pixels>
  requires PixelRange<Pixels>
inline ExpectedMatrixStats ExpectedValues::expected_matrix(const hictk::Chromosome &chrom,
                                                           const hictk::BinTable &bins,
                                                           const Pixels &pixels) const {
  return {pixels,
          chrom,
          chrom,
          bins,
          _expected_weights,
          _expected_scaling_factors.at(chrom),
          *bin_mask(chrom),
          *bin_mask(chrom),
          _min_delta,
          _max_delta};
}

template <typename Pixels>
  requires PixelRange<Pixels>
inline ExpectedMatrixStats ExpectedValues::expected_matrix(const hictk::Chromosome &chrom1,
                                                           const hictk::Chromosome &chrom2,
                                                           const hictk::BinTable &bins,
                                                           const Pixels &pixels) const {
  if (chrom1 == chrom2) {
    return expected_matrix(chrom1, bins, pixels);
  }

  return {pixels,
          chrom1,
          chrom2,
          bins,
          std::vector<N>{},
          1.0,
          *bin_mask(chrom1, chrom2).first,
          *bin_mask(chrom1, chrom2).second,
          _min_delta,
          _max_delta};
}

template <typename File>
  requires HictkSingleResFile<File>
inline auto ExpectedValues::init_pixel_merger_cis(const File &f) {
  using PixelSelector = decltype(std::declval<File>().fetch("chr1"));
  using ThinPixelIt = decltype(std::declval<PixelSelector>().template begin<N>());

  std::vector<PixelSelector> selectors{};

  for (const auto &chrom : f.chromosomes()) {
    if (chrom.is_all()) [[unlikely]] {
      continue;
    }
    auto sel = f.fetch(chrom.name());

    auto first = sel.template begin<N>();
    auto last = sel.template end<N>();

    if (!sel.empty()) [[likely]] {
      selectors.emplace_back(std::move(sel));
    }
  }

  std::vector<ThinPixelIt> heads(selectors.size());
  std::vector<ThinPixelIt> tails(selectors.size());

  for (std::size_t i = 0; i < selectors.size(); ++i) {
    heads[i] = selectors[i].template begin<N>();
    tails[i] = selectors[i].template end<N>();
  }

  return std::make_pair(std::move(selectors),
                        hictk::transformers::PixelMerger{std::move(heads), std::move(tails)});
}

template <typename File>
  requires HictkSingleResFile<File>
inline void ExpectedValues::init_bin_masks(
    const File &f,
    const phmap::flat_hash_map<hictk::Chromosome, std::vector<bool>> &bin_mask_seed) {
  for (const auto &chrom : f.chromosomes()) {
    if (chrom.is_all()) [[unlikely]] {
      continue;
    }
    auto sel = f.fetch(chrom.name());

    auto first = sel.template begin<N>();
    auto last = sel.template end<N>();

    if (!_bin_masks.contains(std::make_pair(chrom, chrom))) {
      if (first == last) {
        SPDLOG_DEBUG("[{}]: no interactions found: masking the entire chromosome!", chrom.name());
        const auto num_bins = (chrom.size() + f.resolution() - 1) / f.resolution();
        add_bin_mask(chrom, std::vector(num_bins, true), bin_mask_seed);
      } else {
        const hictk::transformers::JoinGenomicCoords jsel(first, last, f.bins_ptr());
        add_bin_mask(chrom, mad_max_filtering(jsel, chrom, f.resolution(), _mad_max),
                     bin_mask_seed);
      }
    }
  }
}

}  // namespace nchg
