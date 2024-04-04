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

#include <fmt/format.h>

#include <cstdint>
#include <hictk/file.hpp>
#include <memory>
#include <variant>

#include "nchg/common.hpp"
#include "nchg/config.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {

static void print_header() {
  fmt::print(
      FMT_STRING("chrom1\t"
                 "start1\t"
                 "end2\t"
                 "chrom2\t"
                 "start2\t"
                 "end2\t"
                 "pvalue\t"
                 "observed_count\t"
                 "expected_count\t"
                 "odds_ratio\t"
                 "omega\n"));
}

template <typename FilePtr>
static void process_all_chromosomes(const FilePtr &f, const ComputePvalConfig &c) {
  NCHG nchg(f, c.min_delta, c.max_delta);

  bool header_printed = false;
  for (std::uint32_t chrom1_id = 0; chrom1_id < f->chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = f->chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < f->chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = f->chromosomes().at(chrom2_id);
      nchg.init_matrix(chrom1, chrom2);

      if (c.write_header && !header_printed) {
        print_header();
        header_printed = true;
      }

      nchg.print_pvalues(chrom1, chrom2);
      nchg.erase_matrix(chrom1, chrom2);
    }
  }
}

template <typename FilePtr>
static void process_cis_only_chromosomes(const FilePtr &f, const ComputePvalConfig &c) {
  assert(c.cis_only);
  using File = remove_cvref_t<decltype(*f)>;
  auto nchg = NCHG<File>::cis_only(f, c.min_delta, c.max_delta);

  bool header_printed = false;
  for (const auto &chrom : f->chromosomes()) {
    nchg.init_matrix(chrom);

    if (c.write_header && !header_printed) {
      print_header();
      header_printed = true;
    }

    nchg.print_pvalues(chrom);
    nchg.erase_matrix(chrom);
  }
}

template <typename FilePtr>
static void process_trans_only(const FilePtr &f, const ComputePvalConfig &c) {
  assert(c.trans_only);
  using File = remove_cvref_t<decltype(*f)>;
  auto nchg = NCHG<File>::trans_only(f);

  bool header_printed = false;
  for (std::uint32_t chrom1_id = 0; chrom1_id < f->chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = f->chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < f->chromosomes().size();
         ++chrom2_id) {
      const auto &chrom2 = f->chromosomes().at(chrom2_id);
      nchg.init_matrix(chrom1, chrom2);

      if (c.write_header && !header_printed) {
        print_header();
        header_printed = true;
      }

      nchg.print_pvalues(chrom1, chrom2);
      nchg.erase_matrix(chrom1, chrom2);
    }
  }
}

template <typename FilePtr>
static void process_one_chromosome(const FilePtr &f, const ComputePvalConfig &c) {
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  using File = remove_cvref_t<decltype(*f)>;

  auto nchg = NCHG<File>::chromosome_pair(f, chrom1, chrom2, c.min_delta, c.max_delta);

  nchg.init_matrix(f->chromosomes().at(c.chrom1), f->chromosomes().at(c.chrom2));

  if (c.write_header) {
    print_header();
  }
  nchg.print_pvalues(f->chromosomes().at(c.chrom1), f->chromosomes().at(c.chrom2));
}

int run_nchg_compute(const ComputePvalConfig &c) {
  // clang-format off
  using FilePtr =
      std::variant<std::shared_ptr<const hictk::File>,
          std::shared_ptr<const hictk::cooler::File>,
          std::shared_ptr<const hictk::hic::File>>;
  // clang-format on

  const auto f = [&]() -> FilePtr {
    hictk::File ff(c.path.string(), c.resolution);
    return {std::visit(
        [&](auto &&ff) {
          using FileT = std::remove_reference_t<decltype(ff)>;
          return FilePtr{std::make_shared<const FileT>(std::forward<FileT>(ff))};
        },
        ff.get())};
  }();

  std::visit(
      [&](const auto &f_) {
        if (c.cis_only) {
          assert(c.chrom1 == "all");
          assert(c.chrom2 == "all");
          process_cis_only_chromosomes(f_, c);
          return;
        }

        if (c.trans_only) {
          assert(c.chrom1 == "all");
          assert(c.chrom2 == "all");
          process_trans_only(f_, c);
          return;
        }

        if (c.chrom1 == "all") {
          assert(c.chrom2 == "all");
          process_all_chromosomes(f_, c);
          return;
        }
        process_one_chromosome(f_, c);
      },
      f);

  return 0;
}
}  // namespace nchg
