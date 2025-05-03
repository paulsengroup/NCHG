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

#include "nchg/genomic_domains.hpp"

#include <stdexcept>
#include <utility>

namespace nchg {

BEDPE::BEDPE(const hictk::GenomicInterval& range_) : BEDPE(range_, range_) {}

BEDPE::BEDPE(hictk::GenomicInterval range1_, hictk::GenomicInterval range2_)
    : _range1(std::move(range1_)), _range2(std::move(range2_)) {
  if (_range1 > _range2) {
    throw std::invalid_argument("range1 cannot be greater than range2");
  }
}

bool BEDPE::operator==(const BEDPE& other) const noexcept {
  return _range1 == other._range1 && _range2 == other._range2;
}

bool BEDPE::operator!=(const BEDPE& other) const noexcept { return !(*this == other); }

bool BEDPE::operator<(const BEDPE& other) const noexcept {
  return less_than_op_linear(*this, other);
}

bool BEDPE::less_than_op_linear(const BEDPE& dom1, const BEDPE& dom2) noexcept {
  if (dom1._range1 == dom2._range1) {
    return dom1._range2 < dom2._range2;
  }
  return dom1._range1 < dom2._range1;
}

bool BEDPE::less_than_op_tiled(const BEDPE& dom1, const BEDPE& dom2) noexcept {
  const auto& c1 = dom1._range1.chrom();
  const auto& c2 = dom2._range1.chrom();
  const auto& c3 = dom1._range2.chrom();
  const auto& c4 = dom2._range2.chrom();

  if (c1 != c2) {
    return c1 < c2;
  }
  if (c3 != c4) {
    return c3 < c4;
  }

  auto pos1 = dom1._range1.start();
  auto pos2 = dom2._range1.start();

  if (pos1 != pos2) {
    return pos1 < pos2;
  }

  auto pos3 = dom1._range2.start();
  auto pos4 = dom2._range2.start();

  if (pos3 != pos4) {
    return pos3 < pos4;
  }

  pos1 = dom1._range1.end();
  pos2 = dom2._range1.end();

  if (pos1 != pos2) {
    return pos1 < pos2;
  }

  pos3 = dom1._range2.end();
  pos4 = dom2._range2.end();

  return pos3 < pos4;
}

}  // namespace nchg
