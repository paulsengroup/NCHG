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

#include <cstddef>
#include <functional>
#include <hictk/genomic_interval.hpp>
#include <utility>

#include "nchg/concepts.hpp"
#include "nchg/hash.hpp"

namespace nchg {

constexpr const hictk::GenomicInterval& BEDPE::range1() const noexcept { return _range1; }
constexpr const hictk::GenomicInterval& BEDPE::range2() const noexcept { return _range2; }

template <typename N>
  requires arithmetic<N>
inline BG2Domain<N>::BG2Domain(BEDPE domain_, N count_) noexcept
    : _domain(std::move(domain_)), _count(count_) {}

template <typename N>
  requires arithmetic<N>
constexpr const BEDPE& BG2Domain<N>::domain() const noexcept {
  return _domain;
}

template <typename N>
  requires arithmetic<N>
constexpr N BG2Domain<N>::count() const noexcept {
  return _count;
}

template <typename N>
  requires arithmetic<N>
constexpr auto BG2Domain<N>::operator+=(N n) noexcept -> BG2Domain& {
  _count += n;
  return *this;
}

template <typename N>
  requires arithmetic<N>
inline bool BG2Domain<N>::operator==(const BG2Domain& other) const noexcept {
  return _domain == other._domain && _count == other._count;
}

template <typename N>
  requires arithmetic<N>
inline bool BG2Domain<N>::operator!=(const BG2Domain& other) const noexcept {
  return !(*this == other);
}

template <typename N>
  requires arithmetic<N>
inline bool BG2Domain<N>::operator<(const BG2Domain& other) const noexcept {
  return less_than_op_linear(*this, other);
}

template <typename N>
  requires arithmetic<N>
inline bool BG2Domain<N>::less_than_op_linear(const BG2Domain& dom1,
                                              const BG2Domain& dom2) noexcept {
  if (dom1._domain == dom2._domain) [[unlikely]] {
    return dom1._count < dom2._count;
  }
  return BEDPE::less_than_op_linear(dom1._domain, dom2._domain);
}

template <typename N>
  requires arithmetic<N>
inline bool BG2Domain<N>::less_than_op_tiled(const BG2Domain& dom1,
                                             const BG2Domain& dom2) noexcept {
  if (dom1._domain == dom2._domain) [[unlikely]] {
    return dom1._count < dom2._count;
  }
  return BEDPE::less_than_op_tiled(dom1._domain, dom2._domain);
}

}  // namespace nchg

inline std::size_t std::hash<nchg::BEDPE>::operator()(const nchg::BEDPE& dom) const noexcept {
  return nchg::internal::hash_combine(0, dom.range1(), dom.range2());
}

template <typename N>
inline std::size_t std::hash<nchg::BG2Domain<N>>::operator()(
    const nchg::BG2Domain<N>& dom) const noexcept {
  return nchg::internal::hash_combine(0, dom.domain(), dom.count());
}
