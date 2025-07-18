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

#include <arrow/io/type_fwd.h>
#include <arrow/type_fwd.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <hictk/reference.hpp>
#include <memory>
#include <string>
#include <vector>

namespace parquet {
class StreamReader;
}

namespace nchg {

class ParquetStatsFileReader {
 public:
  enum class RecordType : std::uint_fast8_t { infer, NCHGCompute, NCHGFilter };
  template <typename Stats>
  class iterator;

 private:
  std::filesystem::path _path;
  std::shared_ptr<arrow::Schema> _schema;
  RecordType _type{RecordType::NCHGCompute};
  std::shared_ptr<const hictk::Reference> _chroms;
  std::uint8_t _format_version{1};
  std::string _metadata;
  std::shared_ptr<parquet::StreamReader> _sr;

  ParquetStatsFileReader(std::filesystem::path path,
                         const std::shared_ptr<arrow::io::ReadableFile> &fp, RecordType record_type,
                         std::size_t buffer_size);
  ParquetStatsFileReader(std::filesystem::path path, std::shared_ptr<arrow::io::ReadableFile> fp,
                         std::shared_ptr<const hictk::Reference> chromosomes,
                         RecordType record_type, std::size_t buffer_size);

 public:
  ParquetStatsFileReader() = default;
  ParquetStatsFileReader(const std::filesystem::path &path, RecordType record_type,
                         std::size_t buffer_size = 1'000'000);

  [[nodiscard]] const std::filesystem::path &path() const noexcept;

  [[nodiscard]] std::shared_ptr<arrow::Schema> file_schema() const noexcept;

  [[nodiscard]] auto record_type() const noexcept -> RecordType;

  [[nodiscard]] std::shared_ptr<const hictk::Reference> chromosomes() const noexcept;
  template <typename Stats>
  [[nodiscard]] auto begin() -> iterator<Stats>;
  template <typename Stats>
  [[nodiscard]] auto end() -> iterator<Stats>;

  [[nodiscard]] std::string_view metadata() const noexcept;
  [[nodiscard]] static std::string read_metadata(const std::filesystem::path &path,
                                                 const std::vector<std::string> &ignored_keys = {});
  [[nodiscard]] std::uint8_t format_version() const noexcept;

  template <typename Stats>
  class iterator {
    std::shared_ptr<const hictk::Reference> _chroms;
    std::shared_ptr<parquet::StreamReader> _sr;
    std::shared_ptr<std::string> _buffer;
    Stats _value{};
    std::int64_t _offset{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Stats;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;

    iterator(std::shared_ptr<const hictk::Reference> chroms,
             std::shared_ptr<parquet::StreamReader> sr, bool init_value = true);
    [[nodiscard]] static auto at_end(std::shared_ptr<const hictk::Reference> chroms,
                                     std::shared_ptr<parquet::StreamReader> sr) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const noexcept -> const_reference;
    [[nodiscard]] auto operator->() const noexcept -> const_pointer;

    auto operator++() -> iterator &;

   private:
    void read_pixel();
  };
};

}  // namespace nchg

#include "../../parquet_stats_file_reader_impl.hpp"
