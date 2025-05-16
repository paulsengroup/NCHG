// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
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
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/register/box.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <hictk/chromosome.hpp>
#include <hictk/genomic_interval.hpp>
#include <memory>
#include <numeric>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "nchg/common.hpp"
#include "nchg/concepts.hpp"
#include "nchg/hash.hpp"

namespace nchg::internal {
struct DummyBoostGeometryInterval {
  std::uint32_t start;
  std::uint32_t end;
};
}  // namespace nchg::internal

BOOST_GEOMETRY_REGISTER_POINT_2D(nchg::internal::DummyBoostGeometryInterval, std::uint32_t,
                                 boost::geometry::cs::cartesian, start, end);
BOOST_GEOMETRY_REGISTER_BOX_2D_4VALUES(nchg::BEDPE, nchg::internal::DummyBoostGeometryInterval,
                                       start1(), start2(), end1(), end2());

namespace nchg {

constexpr const hictk::GenomicInterval& BEDPE::range1() const noexcept {
  assert(_range1);
  return *_range1;
}
constexpr const hictk::GenomicInterval& BEDPE::range2() const noexcept {
  assert(_range2);
  return *_range2;
}

constexpr const hictk::Chromosome& BEDPE::chrom1() const noexcept { return range1().chrom(); }

constexpr std::uint32_t BEDPE::start1() const noexcept { return range1().start(); }

constexpr std::uint32_t BEDPE::end1() const noexcept { return range1().end(); }

constexpr const hictk::Chromosome& BEDPE::chrom2() const noexcept { return range2().chrom(); }

constexpr std::uint32_t BEDPE::start2() const noexcept { return range2().start(); }

constexpr std::uint32_t BEDPE::end2() const noexcept { return range2().end(); }

template <typename N>
  requires arithmetic<N>
inline GenomicDomainsIndexed<N>::GenomicDomainsIndexed(hictk::Chromosome chrom1_,
                                                       hictk::Chromosome chrom2_,
                                                       std::span<const BEDPE> domains)
    : GenomicDomainsIndexed(std::move(chrom1_), std::move(chrom2_), build_rtree(domains)) {}

template <typename N>
  requires arithmetic<N>
inline GenomicDomainsIndexed<N>::GenomicDomainsIndexed(hictk::Chromosome chrom1_,
                                                       hictk::Chromosome chrom2_,
                                                       std::shared_ptr<const RTree> domains)
    : _chrom1(std::move(chrom1_)),
      _chrom2(std::move(chrom2_)),
      _domains(std::move(domains)),
      _counts(!!_domains ? _domains->size() : 0, N{0}) {}

template <typename N>
  requires arithmetic<N>
constexpr const hictk::Chromosome& GenomicDomainsIndexed<N>::chrom1() const noexcept {
  return _chrom1;
}

template <typename N>
  requires arithmetic<N>
constexpr const hictk::Chromosome& GenomicDomainsIndexed<N>::chrom2() const noexcept {
  return _chrom2;
}

template <typename N>
  requires arithmetic<N>
inline bool GenomicDomainsIndexed<N>::empty() const noexcept {
  return size() == 0;
}

template <typename N>
  requires arithmetic<N>
inline std::size_t GenomicDomainsIndexed<N>::size() const noexcept {
  if (!_domains) {
    return 0;
  }

  assert(_domains->size() == _counts.size());
  return _domains->size();
}

template <typename N>
  requires arithmetic<N>
inline void GenomicDomainsIndexed<N>::clear_counts() noexcept {
  std::ranges::fill(_counts, N{0});
}

template <typename N>
  requires arithmetic<N>
inline auto GenomicDomainsIndexed<N>::to_vector() const -> std::vector<std::pair<BEDPE, N>> {
  if (empty()) {
    return {};
  }

  std::vector<std::pair<BEDPE, N>> buff(size(), std::make_pair(BEDPE{}, N{}));
  for (const auto& [dom, i] : *_domains) {
    assert(i < size());
    assert(buff[i].first == BEDPE{});
    buff[i] = std::make_pair(BEDPE{{dom.chrom1(), dom.start1(), dom.end1() + 1},
                                   {dom.chrom2(), dom.start2(), dom.end2() + 1}},
                             _counts[i]);
  }

  assert(std::ranges::all_of(buff, [](const auto& dom) { return dom.first != BEDPE{}; }));
  return buff;
}

template <typename N>
  requires arithmetic<N>
inline N GenomicDomainsIndexed<N>::at(const BEDPE& dom) const {
  if (_domains) {
    const BEDPE open_open_dom{{dom.chrom1(), dom.start1(), dom.end1() - 1},
                              {dom.chrom2(), dom.start2(), dom.end2() - 1}};
    auto first = _domains->qbegin(boost::geometry::index::intersects(open_open_dom));
    const auto last = _domains->qend();

    while (first != last) {
      if (first->first == open_open_dom) {
        return _counts[first->second];
      }
      ++first;
    }
  }

  throw std::out_of_range("unable to find domain");
}

template <typename N>
  requires arithmetic<N>
inline N GenomicDomainsIndexed<N>::sum() const noexcept {
  return std::accumulate(_counts.begin(), _counts.end(), N{0});
}

template <typename N>
  requires arithmetic<N>
template <typename M>
inline std::size_t GenomicDomainsIndexed<N>::add_interactions(const hictk::Pixel<M>& p) {
  assert(p.coords.bin1.chrom() == _chrom1);
  assert(p.coords.bin2.chrom() == _chrom2);
  assert(p.count != 0);
  if (empty()) [[unlikely]] {
    return 0;
  }

  const BEDPE query1{{p.coords.bin1.chrom(), p.coords.bin1.start(), p.coords.bin1.end() - 1},
                     {p.coords.bin2.chrom(), p.coords.bin2.start(), p.coords.bin2.end() - 1}};
  const auto count = conditional_static_cast<N>(p.count);

  std::size_t matches = 0;
  auto aggregate_interactions = [&](const auto& domain) {
    ++matches;
    _counts[domain.second] += count;
  };

  std::for_each(_domains->qbegin(boost::geometry::index::intersects(query1)), _domains->qend(),
                aggregate_interactions);

  if (p.coords.bin1.chrom() == p.coords.bin2.chrom() && p.coords.bin1.id() != p.coords.bin2.id()) {
    const BEDPE query2{query1.range2(), query1.range1()};
    std::for_each(_domains->qbegin(boost::geometry::index::intersects(query2)), _domains->qend(),
                  aggregate_interactions);
  }

  return matches;
}

template <typename N>
  requires arithmetic<N>
inline auto GenomicDomainsIndexed<N>::build_rtree(std::span<const BEDPE> domains)
    -> std::shared_ptr<const RTree> {
  if (domains.empty()) {
    return nullptr;
  }

  auto add_index = [i = std::size_t{}](const auto& dom) mutable {
    return std::make_pair(BEDPE{{dom.chrom1(), dom.start1(), dom.end1() - 1},
                                {dom.chrom2(), dom.start2(), dom.end2() - 1}},
                          i++);
  };

  auto indexed_domains = domains | std::views::transform(add_index);
  return std::make_shared<const RTree>(indexed_domains.begin(), indexed_domains.end());
}

template <typename N>
inline auto GenomicDomains::fetch(const hictk::Chromosome& chrom) const
    -> GenomicDomainsIndexed<N> {
  return fetch<N>(chrom, chrom);
}

template <typename N>
inline auto GenomicDomains::fetch(const hictk::Chromosome& chrom1,
                                  const hictk::Chromosome& chrom2) const
    -> GenomicDomainsIndexed<N> {
  return {chrom1, chrom2, equal_range(chrom1, chrom2)};
}

constexpr std::span<BEDPE> GenomicDomains::operator()() noexcept { return _domains; }

constexpr std::span<const BEDPE> GenomicDomains::operator()() const noexcept { return _domains; }

}  // namespace nchg

inline std::size_t std::hash<nchg::BEDPE>::operator()(const nchg::BEDPE& dom) const noexcept {
  return nchg::internal::hash_combine(0, dom.range1(), dom.range2());
}
