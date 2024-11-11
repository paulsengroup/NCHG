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
#include <utility>
#include <vector>

#include "hictk/pixel.hpp"

namespace nchg {

template <typename It>
inline bool KMerger<It>::Node::operator<(const Node &other) const noexcept {
  return value < other.value;
}

template <typename It>
inline bool KMerger<It>::Node::operator>(const Node &other) const noexcept {
  return value > other.value;
}

template <typename It>
inline bool KMerger<It>::Node::operator==(const Node &other) const noexcept {
  return value == other.value;
}

template <typename It>
inline bool KMerger<It>::Node::operator!=(const Node &other) const noexcept {
  return !(*this == other);
}

template <typename It>
inline KMerger<It>::KMerger(const std::vector<It> &heads, const std::vector<It> &tails)
    : KMerger(heads.begin(), heads.end(), tails.begin()) {
  assert(heads.size() == tails.size());
}

template <typename It>
template <typename ItOfIts>
inline KMerger<It>::KMerger(ItOfIts first_head, ItOfIts last_head, ItOfIts first_tail) {
  while (first_head != last_head) {
    if (*first_head != *first_tail) {
      _heads.emplace_back(*first_head);
      _tails.emplace_back(*first_tail);
    }
    ++first_head;
    ++first_tail;
  }
}

template <typename It>
inline auto KMerger<It>::begin() const -> iterator<ItInternal> {
  return iterator<ItInternal>{std::make_shared<std::vector<ItInternal>>(_heads),
                              std::make_shared<std::vector<ItInternal>>(_tails)};
}

template <typename It>
inline auto KMerger<It>::end() const noexcept -> iterator<ItInternal> {
  return {};
}

template <typename It>
inline auto KMerger<It>::read_all() const -> std::vector<T> {
  std::vector<T> buff{};
  std::copy(begin(), end(), std::back_inserter(buff));
  return buff;
}

template <typename It>
template <typename ItInternal>
inline KMerger<It>::iterator<ItInternal>::iterator(const iterator &other)
    : _value(other._value),
      _pqueue(other._pqueue ? std::make_shared<PQueueT>(*other._pqueue) : nullptr),
      _heads(other._heads ? std::make_shared<std::vector<It>>(*other._heads) : nullptr),
      _tails(other._tails ? std::make_shared<std::vector<It>>(*other._tails) : nullptr),
      _i(other._i) {}

template <typename It>
template <typename ItInternal>
inline KMerger<It>::iterator<ItInternal>::iterator(iterator &&other) noexcept
    : _value(other._value),
      _pqueue(std::move(other._pqueue)),
      _heads(std::move(other._heads)),
      _tails(std::move(other._tails)),
      _i(other._i) {}

template <typename It>
template <typename ItInternal>
inline auto KMerger<It>::iterator<ItInternal>::operator=(const iterator &other) -> iterator & {
  if (this == &other) {
    return *this;
  }

  _value = other._value;
  _pqueue = other._pqueue ? std::make_shared<PQueueT>(*other._pqueue) : nullptr;
  _heads = other._heads ? std::make_shared<std::vector<It>>(*other._heads) : nullptr;
  _tails = other._tails ? std::make_shared<std::vector<It>>(*other._tails) : nullptr;
  _i = other._i;

  return *this;
}

template <typename It>
template <typename ItInternal>
inline auto KMerger<It>::iterator<ItInternal>::operator=(iterator &&other) noexcept -> iterator & {
  if (this == &other) {
    return *this;
  }

  _value = std::move(other._value);
  _pqueue = std::move(other._pqueue);
  _heads = std::move(other._heads);
  _tails = std::move(other._tails);
  _i = other._i;

  return *this;
}

template <typename It>
template <typename ItInternal>
inline bool KMerger<It>::iterator<ItInternal>::operator==(const iterator &other) const noexcept {
  if (!_heads || !other._heads) [[unlikely]] {
    // check if we are at end
    return _heads == other._heads;
  }
  return _heads == other._heads && _i == other._i;
}

template <typename It>
template <typename ItInternal>
inline bool KMerger<It>::iterator<ItInternal>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename It>
template <typename ItInternal>
inline auto KMerger<It>::iterator<ItInternal>::operator*() const -> const_reference {
  return *_value;
}

template <typename It>
template <typename ItInternal>
inline auto KMerger<It>::iterator<ItInternal>::operator->() const -> const_pointer {
  return &(**this);
}

template <typename It>
template <typename ItInternal>
inline auto KMerger<It>::iterator<ItInternal>::operator++() -> iterator & {
  assert(_heads);
  _value = next();
  if (!_value) [[unlikely]] {
    _heads = nullptr;
    _tails = nullptr;
  }

  return *this;
}

template <typename It>
template <typename ItInternal>
inline void KMerger<It>::iterator<ItInternal>::replace_top_node() {
  assert(_heads);
  assert(_tails);
  const auto i = _pqueue->top().i;
  _pqueue->pop();
  if (auto &it = (*_heads)[i]; it != (*_tails)[i]) [[likely]] {
    emplace(*it, i);
    std::ignore = ++it;
  }
}

template <typename It>
template <typename ItInternal>
inline void KMerger<It>::iterator<ItInternal>::emplace(T value, std::size_t i) {
#if defined(__apple_build_version__) && __apple_build_version__ < 16000000
  _pqueue->emplace(Node{std::move(value), i});
#else
  _pqueue->emplace(std::move(value), i);
#endif
}

template <typename It>
template <typename ItInternal>
inline auto KMerger<It>::iterator<ItInternal>::next() -> std::optional<T> {
  if (_pqueue->empty()) [[unlikely]] {
    return {};
  }

  auto current_node = _pqueue->top();
  replace_top_node();
  ++_i;
  return current_node.value;
}

}  // namespace nchg
