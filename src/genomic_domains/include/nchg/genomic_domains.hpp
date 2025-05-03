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

#include "nchg/concepts.hpp"

namespace nchg {
class BEDPE {
  hictk::GenomicInterval _range1{};
  hictk::GenomicInterval _range2{};

 public:
  BEDPE() = default;
  explicit BEDPE(const hictk::GenomicInterval &range_);
  BEDPE(hictk::GenomicInterval range1_, hictk::GenomicInterval range2_);

  [[nodiscard]] bool operator==(const BEDPE &other) const noexcept;
  [[nodiscard]] bool operator!=(const BEDPE &other) const noexcept;
  [[nodiscard]] bool operator<(const BEDPE &other) const noexcept;

  [[nodiscard]] static bool less_than_op_tiled(const BEDPE &dom1, const BEDPE &dom2) noexcept;
  [[nodiscard]] static bool less_than_op_linear(const BEDPE &dom1, const BEDPE &dom2) noexcept;

  [[nodiscard]] constexpr const hictk::GenomicInterval &range1() const noexcept;
  [[nodiscard]] constexpr const hictk::GenomicInterval &range2() const noexcept;
};

template <typename N>
  requires arithmetic<N>
class BG2Domain {
  BEDPE _domain{};
  N _count{};

 public:
  BG2Domain() = default;
  explicit BG2Domain(BEDPE domain_, N count_ = 0) noexcept;

  constexpr const BEDPE &domain() const noexcept;
  constexpr N count() const noexcept;

  constexpr auto operator+=(N n) noexcept -> BG2Domain &;

  [[nodiscard]] bool operator==(const BG2Domain &other) const noexcept;
  [[nodiscard]] bool operator!=(const BG2Domain &other) const noexcept;
  [[nodiscard]] bool operator<(const BG2Domain &other) const noexcept;

  [[nodiscard]] static bool less_than_op_tiled(const BG2Domain &dom1,
                                               const BG2Domain &dom2) noexcept;
  [[nodiscard]] static bool less_than_op_linear(const BG2Domain &dom1,
                                                const BG2Domain &dom2) noexcept;
};
}  // namespace nchg

template <>
struct std::hash<nchg::BEDPE> {
  std::size_t operator()(const nchg::BEDPE &dom) const noexcept;
};

template <typename N>
struct std::hash<nchg::BG2Domain<N>> {
  std::size_t operator()(const nchg::BG2Domain<N> &dom) const noexcept;
};

#include "../../genomic_domains_impl.hpp"
