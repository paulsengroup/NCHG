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
// <https://www.gnu.org/licenses/>

#pragma once

#include <cstdint>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <hictk/file.hpp>
#include <hictk/reference.hpp>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "nchg/expected_matrix.hpp"
#include "nchg/genomic_domains.hpp"
#include "nchg/tools/config.hpp"

namespace nchg {

using DomainAggregationStrategy = ComputePvalConfig::DomainAggregationStrategy;

[[nodiscard]] GenomicDomains parse_domains(const hictk::Reference &chroms,
                                           const std::filesystem::path &path, bool keep_cis,
                                           bool keep_trans,
                                           const std::optional<hictk::Chromosome> &chrom1 = {},
                                           const std::optional<hictk::Chromosome> &chrom2 = {});

[[nodiscard]] GenomicDomains parse_domains(const hictk::Reference &chroms,
                                           const std::filesystem::path &path, bool keep_cis,
                                           bool keep_trans,
                                           const std::optional<std::string> &chrom1 = {},
                                           const std::optional<std::string> &chrom2 = {});
[[nodiscard]] std::filesystem::path write_domains_to_file(const GenomicDomains &domains,
                                                          const std::filesystem::path &dest_dir,
                                                          const hictk::Chromosome &chrom1,
                                                          const hictk::Chromosome &chrom2,
                                                          bool force);

void write_chrom_sizes_to_file(const hictk::Reference &chroms, const std::filesystem::path &path,
                               bool force);

[[nodiscard]] std::vector<std::tuple<BEDPE, std::uint64_t, double>> map_interactions_to_domains(
    const hictk::File &f, const GenomicDomains &domains, const ExpectedMatrixStats &expected_matrix,
    const hictk::Chromosome &chrom1, const hictk::Chromosome &chrom2, std::uint64_t min_delta,
    std::uint64_t max_delta, const std::vector<bool> &bin1_mask, const std::vector<bool> &bin2_mask,
    DomainAggregationStrategy aggregation_stategy);
}  // namespace nchg
