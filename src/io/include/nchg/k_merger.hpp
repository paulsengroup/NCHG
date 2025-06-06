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

#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <queue>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace nchg {

/// This class is basically a wrapper around a priority queue of objects of type Node, and can be
/// used to merge two or more sorted sequences of objects.
/// Nodes consist of an object and an index.
/// The object should be comparable with objects of the same type to determine which comes first in
/// a sorted sequence.
/// The index represents from which iterator the object was read.
/// This allows us to know from which iterator we should read the next object (i.e. the same
/// iterator from which the top node originated)
template <typename It>
class KMerger {
  using ItInternal = std::remove_cvref_t<It>;
  using T = std::remove_cvref_t<decltype(*std::declval<It>())>;
  struct Node {
    T value{};        // NOLINT
    std::size_t i{};  // NOLINT

    bool operator<(const Node &other) const noexcept;
    bool operator>(const Node &other) const noexcept;
    bool operator==(const Node &other) const noexcept;
    bool operator!=(const Node &other) const noexcept;
  };

  std::vector<ItInternal> _heads{};
  std::vector<ItInternal> _tails{};

 public:
  template <typename ItInternal = ItInternal>
  class iterator;

  KMerger() = delete;
  KMerger(const std::vector<It> &heads, const std::vector<It> &tails);

  template <typename ItOfIts>
  KMerger(ItOfIts first_head, ItOfIts last_head, ItOfIts first_tail);

  [[nodiscard]] auto begin() const -> iterator<ItInternal>;
  [[nodiscard]] auto end() const noexcept -> iterator<ItInternal>;

  [[nodiscard]] auto read_all() const -> std::vector<T>;

  template <typename ItInternal>
  class iterator {
    std::optional<T> _value{};

    using PQueueT = std::priority_queue<Node, std::vector<Node>, std::greater<>>;
    std::shared_ptr<PQueueT> _pqueue{};

    std::shared_ptr<std::vector<ItInternal>> _heads{};
    std::shared_ptr<std::vector<ItInternal>> _tails{};
    std::size_t _i{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::input_iterator_tag;

    iterator() = default;

    // apple-clang < 16 chokes if this ctor is implemented out of line
    iterator(std::shared_ptr<std::vector<ItInternal>> heads,
             std::shared_ptr<std::vector<ItInternal>> tails)
        : _pqueue(std::make_shared<PQueueT>()), _heads(std::move(heads)), _tails(std::move(tails)) {
      assert(_heads->size() == _tails->size());
      for (auto &it : *_heads) {
        emplace(*it, _pqueue->size());
        std::ignore = ++it;
      }
      _value = next();
    }

    iterator(const iterator &other);
    iterator(iterator &&other) noexcept;
    ~iterator() noexcept = default;

    auto operator=(const iterator &other) -> iterator &;
    auto operator=(iterator &&other) noexcept -> iterator &;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    auto operator*() const -> const_reference;
    auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;

   private:
    [[nodiscard]] auto next() -> std::optional<T>;
    void replace_top_node();
    void emplace(T value, std::size_t i);
  };
};
}  // namespace nchg

#include "../k_merger_impl.hpp"  // NOLINT
