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
#include <fmt/std.h>

#include <cstdint>
#include <fstream>
#include <hictk/file.hpp>
#include <hictk/fmt/pixel.hpp>
#include <hictk/genomic_interval.hpp>
#include <hictk/reference.hpp>
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

[[nodiscard]] static std::string_view truncate_bed3_record(std::string_view record,
                                                           char sep = '\t') {
  const auto pos1 = record.find('\t');
  if (pos1 == std::string_view::npos) {
    throw std::runtime_error("invalid bed record, expected 3 tokens, found 1");
  }
  const auto pos2 = record.find('\t', pos1 + 1);
  if (pos2 == std::string_view::npos) {
    throw std::runtime_error("invalid bed record, expected 3 tokens, found 2");
  }
  const auto pos3 = record.find('\t', pos2 + 1);

  return record.substr(pos3);
}

[[nodiscard]] static std::vector<hictk::GenomicInterval> parse_domains(
    const hictk::Reference &chroms, const std::filesystem::path &path, std::string_view chrom1,
    std::string_view chrom2) {
  std::vector<hictk::GenomicInterval> domains{};
  std::string buffer{};

  std::ifstream fs{};
  fs.exceptions(fs.exceptions() | std::ios::badbit | std::ios::failbit);

  try {
    fs.open(path);

    for (std::size_t i = 1; std::getline(fs, buffer); ++i) {
      if (buffer.empty()) {
        continue;
      }

      try {
        const auto record = truncate_bed3_record(buffer);
        auto domain = hictk::GenomicInterval::parse_bed(chroms, record);

        if (chrom1 != "all") {
          assert(chrom2 != "all");
          if (domain.chrom().name() != chrom1 && domain.chrom().name() != chrom2) {
            continue;
          }
        }

        domains.emplace_back(std::move(domain));
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("found an invalid record at line {} of file {}: {}"), i, path, e.what()));
      }
    }

  } catch (const std::exception &e) {
    if (!fs.eof()) {
      throw;
    }
  }

  std::sort(domains.begin(), domains.end());
  return domains;
}

template <typename FilePtr, typename File = remove_cvref_t<decltype(*std::declval<FilePtr>())>>
[[nodiscard]] static NCHG<File> init_nchg(const FilePtr &f, const ComputePvalConfig &c) {
  if (c.cis_only) {
    return NCHG<File>::cis_only(f, c.min_delta, c.max_delta);
  }
  if (c.trans_only) {
    return NCHG<File>::trans_only(f);
  }

  NCHG nchg(f, c.min_delta, c.max_delta);
  if (c.chrom1 != "all") {
    assert(c.chrom2 != "all");
    nchg.init_matrix(f->chromosomes().at(c.chrom1), f->chromosomes().at(c.chrom2));
  } else {
    nchg.init_matrices();
  }
  return nchg;
}

template <typename FilePtr>
static void process_domains(const FilePtr &f, const ComputePvalConfig &c) {
  assert(std::filesystem::exists(c.path_to_domains));

  const auto domains = parse_domains(f->chromosomes(), c.path_to_domains, c.chrom1, c.chrom2);

  if (domains.empty()) {
    return;
  }

  const auto nchg = init_nchg(f, c);

  if (c.write_header) {
    print_header();
  }

  for (std::size_t i = 0; i < domains.size(); ++i) {
    for (std::size_t j = i; j < domains.size(); ++j) {
      const auto &d1 = domains[i];
      const auto &d2 = domains[j];

      if (c.cis_only && d1.chrom() != d2.chrom()) {
        break;
      }

      if (c.trans_only && d1.chrom() == d2.chrom()) {
        continue;
      }

      const auto s = nchg.compute(d1, d2);
      fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval, s.pixel.count,
                 s.expected, s.odds_ratio, s.omega);
    }
  }
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

      std::for_each(nchg.begin(chrom1, chrom2), nchg.end(chrom1, chrom2), [&](const auto &s) {
        fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval,
                   s.pixel.count, s.expected, s.odds_ratio, s.omega);
      });
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

    std::for_each(nchg.begin(chrom, chrom), nchg.end(chrom, chrom), [&](const auto &s) {
      fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval, s.pixel.count,
                 s.expected, s.odds_ratio, s.omega);
    });
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

      std::for_each(nchg.begin(chrom1, chrom2), nchg.end(chrom1, chrom2), [&](const auto &s) {
        fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval,
                   s.pixel.count, s.expected, s.odds_ratio, s.omega);
      });
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

  nchg.init_matrix(chrom1, chrom2);

  if (c.write_header) {
    print_header();
  }

  std::for_each(nchg.begin(chrom1, chrom2), nchg.end(chrom1, chrom2), [&](const auto &s) {
    fmt::print(FMT_COMPILE("{:bg2}\t{}\t{}\t{}\t{}\t{}\n"), s.pixel.coords, s.pval, s.pixel.count,
               s.expected, s.odds_ratio, s.omega);
  });
}

int run_nchg_compute(const ComputePvalConfig &c) {
  // clang-format off
  using FilePtr =
      std::variant<
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
        if (!c.path_to_domains.empty()) {
          process_domains(f_, c);
          return;
        }

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
