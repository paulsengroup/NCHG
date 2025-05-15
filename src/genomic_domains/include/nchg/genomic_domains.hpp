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

#include <boost/geometry/index/rtree.hpp>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <hictk/chromosome.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/pixel.hpp>
#include <hictk/reference.hpp>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "nchg/concepts.hpp"

namespace nchg {
class BEDPE {
  std::shared_ptr<const hictk::GenomicInterval> _range1{};
  std::shared_ptr<const hictk::GenomicInterval> _range2{};

  inline static const auto _null_range = std::make_shared<const hictk::GenomicInterval>();

 public:
  BEDPE();
  explicit BEDPE(hictk::GenomicInterval range_);
  explicit BEDPE(const std::shared_ptr<const hictk::GenomicInterval> &range_);
  BEDPE(hictk::GenomicInterval range1_, hictk::GenomicInterval range2_);
  BEDPE(std::shared_ptr<const hictk::GenomicInterval> range1_,
        std::shared_ptr<const hictk::GenomicInterval> range2_);

  [[nodiscard]] bool operator==(const BEDPE &other) const noexcept;
  [[nodiscard]] bool operator!=(const BEDPE &other) const noexcept;
  [[nodiscard]] bool operator<(const BEDPE &other) const noexcept;

  [[nodiscard]] static bool less_than_op_tiled(const BEDPE &dom1, const BEDPE &dom2) noexcept;
  [[nodiscard]] static bool less_than_op_linear(const BEDPE &dom1, const BEDPE &dom2) noexcept;

  [[nodiscard]] constexpr const hictk::GenomicInterval &range1() const noexcept;
  [[nodiscard]] constexpr const hictk::GenomicInterval &range2() const noexcept;

  [[nodiscard]] constexpr const hictk::Chromosome &chrom1() const noexcept;
  [[nodiscard]] constexpr std::uint32_t start1() const noexcept;
  [[nodiscard]] constexpr std::uint32_t end1() const noexcept;
  [[nodiscard]] constexpr const hictk::Chromosome &chrom2() const noexcept;
  [[nodiscard]] constexpr std::uint32_t start2() const noexcept;
  [[nodiscard]] constexpr std::uint32_t end2() const noexcept;

  // These accessors are required in order to register BEDPE with boost::geometry.
  // They do nothing and throw an exception when called.
  [[noreturn]] constexpr hictk::Chromosome &chrom1();
  [[noreturn]] constexpr std::uint32_t &start1();
  [[noreturn]] constexpr std::uint32_t &end1();
  [[noreturn]] constexpr hictk::Chromosome &chrom2();
  [[noreturn]] constexpr std::uint32_t &start2();
  [[noreturn]] constexpr std::uint32_t &end2();
};

template <typename N>
  requires arithmetic<N>
class GenomicDomainsIndexed {
  hictk::Chromosome _chrom1{};
  hictk::Chromosome _chrom2{};

  using Value = std::pair<BEDPE, std::size_t>;
  using RTree = boost::geometry::index::rtree<Value, boost::geometry::index::rstar<16>>;

  std::shared_ptr<const RTree> _domains{};
  std::vector<N> _counts{};

 public:
  GenomicDomainsIndexed() = default;
  GenomicDomainsIndexed(hictk::Chromosome chrom1_, hictk::Chromosome chrom2_,
                        std::span<const BEDPE> domains);
  GenomicDomainsIndexed(hictk::Chromosome chrom1_, hictk::Chromosome chrom2_,
                        std::shared_ptr<const RTree> domains);

  [[nodiscard]] constexpr const hictk::Chromosome &chrom1() const noexcept;
  [[nodiscard]] constexpr const hictk::Chromosome &chrom2() const noexcept;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void clear_counts() noexcept;

  [[nodiscard]] auto to_vector() const -> std::vector<std::pair<BEDPE, N>>;
  [[nodiscard]] N at(const BEDPE &dom) const;
  [[nodiscard]] N sum() const noexcept;

  template <typename M>
  std::size_t add_interactions(const hictk::Pixel<M> &p);

 private:
  [[nodiscard]] static auto build_rtree(std::span<const BEDPE> domains)
      -> std::shared_ptr<const RTree>;
};

class GenomicDomains {
  std::vector<BEDPE> _domains;

 public:
  GenomicDomains() = default;
  explicit GenomicDomains(std::vector<BEDPE> domains_, bool sort = true);

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool contains(const hictk::Chromosome &chrom) const noexcept;
  [[nodiscard]] bool contains(const hictk::Chromosome &chrom1,
                              const hictk::Chromosome &chrom2) const noexcept;

  template <typename N>
  [[nodiscard]] auto fetch(const hictk::Chromosome &chrom) const -> GenomicDomainsIndexed<N>;
  template <typename N>
  [[nodiscard]] auto fetch(const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2) const
      -> GenomicDomainsIndexed<N>;

  [[nodiscard]] constexpr std::span<BEDPE> operator()() noexcept;
  [[nodiscard]] constexpr std::span<const BEDPE> operator()() const noexcept;

 private:
  [[nodiscard]] std::span<const BEDPE> equal_range(const hictk::Chromosome &chrom) const noexcept;
  [[nodiscard]] std::span<const BEDPE> equal_range(const hictk::Chromosome &chrom1,
                                                   const hictk::Chromosome &chrom2) const noexcept;
};
}  // namespace nchg

template <>
struct std::hash<nchg::BEDPE> {
  std::size_t operator()(const nchg::BEDPE &dom) const noexcept;
};

#include "../../genomic_domains_impl.hpp"
