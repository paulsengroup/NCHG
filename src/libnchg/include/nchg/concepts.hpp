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

#include <concepts>
#include <hictk/file.hpp>
#include <hictk/pixel.hpp>
#include <iterator>
#include <type_traits>

#include "nchg/type_traits.hpp"

namespace nchg {

template <typename T>
concept arithmetic = std::integral<T> || std::floating_point<T>;

template <typename It>
concept PixelStream = requires(It it) {
  std::input_iterator<It>;
  { it->coords } -> std::convertible_to<hictk::PixelCoordinates>;
  { std::remove_cvref_t<decltype(it->count)>{} } -> arithmetic;
};

template <typename It>
concept ThinPixelStream = requires(It it) {
  std::input_iterator<It>;
  { it->bin1_id } -> std::integral;
  { it->bin2_id } -> std::integral;
  { it->count } -> arithmetic;
};

template <typename File>
concept HictkSingleResFile =
    // clang-format off
    std::is_same_v<hictk::File, File>      ||
    std::is_same_v<hictk::hic::File, File> ||
    std::is_same_v<hictk::cooler::File, File>
    // clang-format on
    ;

template <typename T>
concept UniquePtr = is_unique_ptr_v<T>;
template <typename T>
concept SharedPtr = is_shared_ptr_v<T>;
template <typename T>
concept SmartPtr = UniquePtr<T> || SharedPtr<T>;

}  // namespace nchg
