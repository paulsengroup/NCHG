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

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <hictk/genomic_interval.hpp>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace nchg {

BEDPE::BEDPE() : BEDPE(_null_range) {}

BEDPE::BEDPE(hictk::GenomicInterval range_)
    : BEDPE(std::make_shared<const hictk::GenomicInterval>(std::move(range_))) {}

BEDPE::BEDPE(const std::shared_ptr<const hictk::GenomicInterval>& range_) : BEDPE(range_, range_) {}

BEDPE::BEDPE(hictk::GenomicInterval range1_, hictk::GenomicInterval range2_)
    : BEDPE(std::make_shared<const hictk::GenomicInterval>(std::move(range1_)),
            std::make_shared<const hictk::GenomicInterval>(std::move(range2_))) {}

BEDPE::BEDPE(std::shared_ptr<const hictk::GenomicInterval> range1_,
             std::shared_ptr<const hictk::GenomicInterval> range2_)
    : _range1(std::move(range1_)), _range2(std::move(range2_)) {
  assert(_range1);
  assert(_range2);
}

bool BEDPE::operator==(const BEDPE& other) const noexcept {
  assert(_range1);
  assert(_range2);
  return *_range1 == *other._range1 && *_range2 == *other._range2;
}

bool BEDPE::operator!=(const BEDPE& other) const noexcept { return !(*this == other); }

bool BEDPE::operator<(const BEDPE& other) const noexcept {
  return less_than_op_linear(*this, other);
}

bool BEDPE::less_than_op_linear(const BEDPE& dom1, const BEDPE& dom2) noexcept {
  if (dom1.range1() == dom2.range1()) {
    return dom1.range2() < dom2.range2();
  }
  return dom1.range1() < dom2.range1();
}

bool BEDPE::less_than_op_tiled(const BEDPE& dom1, const BEDPE& dom2) noexcept {
  const auto& c1 = dom1.range1().chrom();
  const auto& c2 = dom2.range1().chrom();
  const auto& c3 = dom1.range2().chrom();
  const auto& c4 = dom2.range2().chrom();

  if (c1 != c2) {
    return c1 < c2;
  }
  if (c3 != c4) {
    return c3 < c4;
  }

  auto pos1 = dom1.range1().start();
  auto pos2 = dom2.range1().start();

  if (pos1 != pos2) {
    return pos1 < pos2;
  }

  auto pos3 = dom1.range2().start();
  auto pos4 = dom2.range2().start();

  if (pos3 != pos4) {
    return pos3 < pos4;
  }

  pos1 = dom1.range1().end();
  pos2 = dom2.range1().end();

  if (pos1 != pos2) {
    return pos1 < pos2;
  }

  pos3 = dom1.range2().end();
  pos4 = dom2.range2().end();

  return pos3 < pos4;
}

GenomicDomains::GenomicDomains(std::vector<BEDPE> domains_, bool sort)
    : _domains(std::move(domains_)) {
  if (sort) {
    std::ranges::sort(_domains, BEDPE::less_than_op_tiled);
  } else {
    assert(std::ranges::is_sorted(_domains, BEDPE::less_than_op_tiled));
  }

  const auto [first, last] = std::ranges::unique(_domains);
  const auto num_duplicates = std::distance(first, last);
  if (num_duplicates != 0) {
    _domains.erase(first, last);
    SPDLOG_WARN("dropped {} duplicated domain(s)", num_duplicates);
  }
}

bool GenomicDomains::empty() const noexcept { return size() == 0; }

std::size_t GenomicDomains::size() const noexcept { return _domains.size(); }

bool GenomicDomains::contains(const hictk::Chromosome& chrom) const noexcept {
  return contains(chrom, chrom);
}

bool GenomicDomains::contains(const hictk::Chromosome& chrom1,
                              const hictk::Chromosome& chrom2) const noexcept {
  return !equal_range(chrom1, chrom2).empty();
}

std::span<const BEDPE> GenomicDomains::equal_range(const hictk::Chromosome& chrom) const noexcept {
  return equal_range(chrom, chrom);
}

std::span<const BEDPE> GenomicDomains::equal_range(const hictk::Chromosome& chrom1,
                                                   const hictk::Chromosome& chrom2) const noexcept {
  auto cmp = [](const auto& dom1, const auto& dom2) noexcept {
    if (dom1.range1().chrom() == dom2.range1().chrom()) {
      return dom1.range2().chrom() < dom2.range2().chrom();
    }
    return dom1.range1().chrom() < dom2.range1().chrom();
  };

  const BEDPE query{hictk::GenomicInterval{chrom1}, hictk::GenomicInterval{chrom2}};
  const auto [first, last] = std::ranges::equal_range(_domains, query, cmp);
  return {first, last};
}

// NOLINTBEGIN(*-member-functions-to-static)
hictk::Chromosome& BEDPE::chrom1() { throw std::runtime_error("not implemented!"); }

std::uint32_t& BEDPE::start1() { chrom1(); }

std::uint32_t& BEDPE::end1() { start1(); }

hictk::Chromosome& BEDPE::chrom2() { chrom1(); }

std::uint32_t& BEDPE::start2() { start1(); }

std::uint32_t& BEDPE::end2() { end1(); }
// NOLINTEND(*-member-functions-to-static)

}  // namespace nchg
