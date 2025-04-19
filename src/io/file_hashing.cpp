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

// clang-format off
#include "nchg/suppress_warnings.hpp"
// clang-format on

#include "nchg/file_hashing.hpp"

#include <fmt/format.h>
NCHG_DISABLE_WARNING_PUSH
NCHG_DISABLE_WARNING_MACRO_REDEFINED
NCHG_DISABLE_WARNING_OLD_STYLE_CAST
#include <xxh3.h>
NCHG_DISABLE_WARNING_POP

#include <cassert>
#include <exception>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

namespace nchg {

XXH3Digest::XXH3Digest() : _digest(fmt::format("{:0>{}}", "", WIDTH / 4)) {}
XXH3Digest::XXH3Digest(std::string digest) : _digest(std::move(digest)) {
  static_assert(WIDTH != 0);
  static_assert(WIDTH % 4 == 0);
  if (_digest.size() != WIDTH / 4) {
    throw std::invalid_argument(
        fmt::format("digest should be exactly {} characters long", WIDTH / 4));
  }

  constexpr std::string_view valid_hex_chars = "0123456789abcdefABCDEF";
  if (_digest.find_first_not_of(valid_hex_chars) != std::string::npos) {
    throw std::invalid_argument("digest is not a valid hex string");
  }
}

std::string_view XXH3Digest::operator()() const noexcept { return _digest; }

void XXH3StateDeleter::operator()(XXH3_state_t* state) const noexcept {
  if (state) {
    XXH3_freeState(state);
  }
}

XXH3StatePtr init_xxh3_state() {
  XXH3StatePtr state(XXH3_createState());
  if (!state) {
    throw std::bad_alloc();
  }
  const auto status = XXH3_128bits_reset(state.get());

  if (status == XXH_ERROR) {
    throw std::runtime_error("hashing failed: XXH3_128bits_reset() failed");
  }
  return state;
}

void update_state(XXH3StatePtr& state, std::string_view buffer) {
  const auto status = XXH3_128bits_update(state.get(), buffer.data(), buffer.size());
  if (status == XXH_ERROR) [[unlikely]] {
    throw std::runtime_error("hashing failed: XXH3_128bits_update() failed");
  }
}

[[nodiscard]] std::string to_hex(const XXH3StatePtr& state) {
  assert(state);
  const auto digest = XXH3_128bits_digest(state.get());
  return fmt::format("{:016x}{:016x}", digest.high64, digest.low64);
}

[[nodiscard]] static XXH3StatePtr digest_chunk(std::istream& fs, const std::streamoff& pos,
                                               const std::streamsize& size,
                                               XXH3StatePtr state = nullptr,
                                               // NOLINTNEXTLINE
                                               std::size_t chunk_size = 4ULL << 20ULL) {
  if (!state) {
    state = init_xxh3_state();
  }

  try {
    fs.seekg(pos, std::ios::beg);
    std::string buffer{};

    const auto end_pos = pos + size;
    while (fs.tellg() < end_pos) {
      buffer.resize(std::min(chunk_size, static_cast<std::size_t>(end_pos - fs.tellg())));
      fs.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      update_state(state, buffer);
    }
  } catch (...) {
    if (!fs.eof()) {
      throw;
    }
  }

  return state;
}

std::string hash_file(const std::filesystem::path& path, std::streamsize sample_size) {
  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  try {
    if (sample_size < 1) {
      throw std::invalid_argument("sample_size must be greater than zero");
    }

    fs.open(path, std::ios::in | std::ios::binary | std::ios::ate);

    const auto file_size = static_cast<std::streamsize>(fs.tellg());
    const auto head_size =
        sample_size >= file_size ? file_size : std::max(sample_size / 2, std::streamsize{1});
    const auto tail_size = sample_size >= file_size ? 0 : sample_size - head_size;

    auto state = digest_chunk(fs, 0, head_size);
    if (tail_size != 0) {
      assert(state);
      const auto pos = file_size - tail_size;
      state = digest_chunk(fs, pos, tail_size, std::move(state));
    }

    return to_hex(state);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format("failed to hash file \"{}\": {}", path.string(), e.what()));
  } catch (...) {
    throw std::runtime_error(
        fmt::format("failed to hash file \"{}\": unknown error", path.string()));
  }
}

}  // namespace nchg
