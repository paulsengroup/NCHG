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

#include <cstdint>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/transformers/join_genomic_coords.hpp>
#include <highfive/H5File.hpp>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/nchg.hpp"
#include "nchg/tools.hpp"

namespace nchg {

template <typename FilePtr>
[[nodiscard]] static std::pair<std::vector<double>, phmap::btree_map<hictk::Chromosome, double>>
compute_cis_profiles(FilePtr f, std::uint64_t min_delta, std::uint64_t max_delta,
                     std::uint64_t num_quantiles) {
  return NCHG{f, min_delta, max_delta, num_quantiles}.compute_expected_profile();
}

template <typename FilePtr>
[[nodiscard]] static phmap::btree_map<std::pair<hictk::Chromosome, hictk::Chromosome>, double>
compute_nnz_avg_values(FilePtr f) {
  phmap::btree_map<std::pair<hictk::Chromosome, hictk::Chromosome>, double> weights{};
  for (std::uint32_t chrom1_id = 0; chrom1_id < f->chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = f->chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < f->chromosomes().size();
         ++chrom2_id) {
      const auto &chrom2 = f->chromosomes().at(chrom2_id);

      SPDLOG_INFO(FMT_STRING("processing {}:{}..."), chrom1.name(), chrom2.name());

      const auto sel = f->fetch(chrom1.name(), chrom2.name());
      const hictk::transformers::JoinGenomicCoords jsel(
          sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(), f->bins_ptr());

      weights.emplace(
          std::make_pair(chrom1, chrom2),
          ExpectedMatrix{jsel.begin(), jsel.end(), chrom1, chrom2, f->bins(), std::vector<double>{},
                         0, std::numeric_limits<std::uint64_t>::max()}

              .nnz_avg());
    }
  }

  return weights;
}

static void write_chromosomes(HighFive::File &f, const hictk::Reference &chromosomes) {
  auto grp = f.createGroup("chroms");

  std::vector<std::string> chromosome_names(chromosomes.size());
  std::transform(chromosomes.begin(), chromosomes.end(), chromosome_names.begin(),
                 [](const hictk::Chromosome &c) { return std::string{c.name()}; });
  grp.createDataSet("name", chromosome_names);

  std::vector<std::uint32_t> chromosome_sizes(chromosomes.size());
  std::transform(chromosomes.begin(), chromosomes.end(), chromosome_sizes.begin(),
                 [](const hictk::Chromosome &c) { return c.size(); });
  grp.createDataSet("length", chromosome_sizes);
}

static void write_profile(HighFive::File &f, const std::vector<double> &profile,
                          const phmap::btree_map<hictk::Chromosome, double> &scaling_factors) {
  auto grp = f.createGroup("profile");
  if (!profile.empty()) {
    HighFive::DataSetCreateProps props{};
    props.add(HighFive::Chunking({profile.size()}));
    props.add(HighFive::Deflate(9));
    grp.createDataSet("values", profile, props);

    std::vector<double> scaling_factors_flat(scaling_factors.size());
    std::transform(scaling_factors.begin(), scaling_factors.end(), scaling_factors_flat.begin(),
                   [](const auto &kv) { return kv.second; });
    grp.createDataSet("scaling-factors", scaling_factors_flat);
  }
}

static void write_nnz_avg_values(
    HighFive::File &f,
    const phmap::btree_map<std::pair<hictk::Chromosome, hictk::Chromosome>, double>
        &nnz_avg_values) {
  auto grp = f.createGroup("avg-values");
  if (nnz_avg_values.empty()) {
    return;
  }

  std::vector<std::string> chrom1{};
  std::vector<std::string> chrom2{};
  std::vector<double> values{};

  for (const auto &[cp, value] : nnz_avg_values) {
    chrom1.emplace_back(std::string{cp.first.name()});   // NOLINT
    chrom2.emplace_back(std::string{cp.second.name()});  // NOLINT
    values.emplace_back(value);
  }

  grp.createDataSet("chrom1", chrom1);
  grp.createDataSet("chrom2", chrom2);
  grp.createDataSet("value", values);
}

static void write_expected_values(
    const std::filesystem::path &path, std::string_view input_file_basename,
    std::uint32_t resolution, const hictk::Reference &chromosomes,
    const std::vector<double> &profile,
    const phmap::btree_map<hictk::Chromosome, double> &scaling_factors,
    const phmap::btree_map<std::pair<hictk::Chromosome, hictk::Chromosome>, double>
        &nnz_avg_values) {
  SPDLOG_INFO(FMT_STRING("writing expected value profiles to {}..."), path);

  HighFive::File f(path.string(), HighFive::File::Create);
  f.createAttribute("source-file", std::string{input_file_basename});
  f.createAttribute("resolution", resolution);

  write_chromosomes(f, chromosomes);
  write_profile(f, profile, scaling_factors);
  write_nnz_avg_values(f, nnz_avg_values);
}

template <typename FilePtr>
static void process_all_chromosomes(FilePtr f, const ExpectedConfig &c) {
  const auto &[profile, scaling_factors] =
      compute_cis_profiles(f, c.min_delta, c.max_delta, c.num_quantiles);
  const auto nnz_avg_values = compute_nnz_avg_values(f);

  if (c.force) {
    std::filesystem::remove(c.output_path);
  }

  write_expected_values(c.output_path, std::filesystem::path{f->path()}.filename().string(),
                        f->resolution(), f->chromosomes(), profile, scaling_factors,
                        nnz_avg_values);
}

template <typename FilePtr>
static void process_one_chromosome_pair(FilePtr f, const ExpectedConfig &c) {
  assert(c.chrom1 != "all");
  const auto &chrom1 = f->chromosomes().at(c.chrom1);
  const auto &chrom2 = f->chromosomes().at(c.chrom2);

  SPDLOG_INFO(FMT_STRING("processing {}:{}..."), chrom1.name(), chrom2.name());

  const auto sel1 = f->fetch();
  const auto sel2 = f->fetch(chrom1.name(), chrom2.name());

  const hictk::transformers::JoinGenomicCoords jsel1(
      sel1.template begin<std::uint32_t>(), sel1.template end<std::uint32_t>(), f->bins_ptr());
  const hictk::transformers::JoinGenomicCoords jsel2(
      sel2.template begin<std::uint32_t>(), sel2.template end<std::uint32_t>(), f->bins_ptr());

  using PixelIt = decltype(jsel2.begin());
  const ExpectedMatrix<PixelIt> m(jsel2.begin(), jsel2.end(), jsel1.begin(), jsel1.end(), chrom1,
                                  chrom2, f->bins(), c.min_delta, c.max_delta);

  phmap::btree_map<std::pair<hictk::Chromosome, hictk::Chromosome>, double> nnz_avg_values{};
  phmap::btree_map<hictk::Chromosome, double> scaling_factors;
  if (chrom1 != chrom2) {
    nnz_avg_values.emplace(std::make_pair(chrom1, chrom2), m.nnz_avg());
  } else {
    scaling_factors.emplace(chrom1, 1.0);
  }

  if (c.force) {
    std::filesystem::remove(c.output_path);
  }

  write_expected_values(c.output_path, std::filesystem::path{f->path()}.filename().string(),
                        f->resolution(), f->chromosomes(), m.weights(), scaling_factors,
                        nnz_avg_values);
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
