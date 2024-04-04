// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#include <parallel_hashmap/btree.h>
#include "nchg/common.hpp"

#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <highfive/H5File.hpp>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/expected_values.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {

template <typename FilePtr>
static void process_all_chromosomes(FilePtr f, const ExpectedConfig &c) {
  const ExpectedValues evs(f);
  if (c.force) {
    std::filesystem::remove(c.output_path);
  }
  evs.serialize(c.output_path);
}

template <typename FilePtr>
static void process_one_chromosome_pair(FilePtr f, const ExpectedConfig &c) {
  assert(c.chrom1 != "all");
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  using File = remove_cvref_t<decltype(*std::declval<FilePtr>())>;
  const auto evs = ExpectedValues<File>::chromosome_pair(f, chrom1, chrom2);
  if (c.force) {
    std::filesystem::remove(c.output_path);
  }
  evs.serialize(c.output_path);
}

int run_nchg_expected(const ExpectedConfig &c) {
  // clang-format off
  using FilePtr =
      std::variant<std::shared_ptr<const hictk::File>,
          std::shared_ptr<const hictk::cooler::File>,
          std::shared_ptr<const hictk::hic::File>>;
  // clang-format on

  const auto f = [&]() -> FilePtr {
    hictk::File ff(c.input_path.string(), c.resolution);
    return {std::visit(
        [&](auto &&ff) {
          using FileT = std::remove_reference_t<decltype(ff)>;
          return FilePtr{std::make_shared<const FileT>(std::forward<FileT>(ff))};
        },
        ff.get())};
  }();

  std::visit(
      [&](const auto &f_) {
        if (c.chrom1 == "all") {
          assert(c.chrom2 == "all");
          process_all_chromosomes(f_, c);
          return;
        }
        process_one_chromosome_pair(f_, c);
      },
      f);

  return 0;
}

}  // namespace nchg
