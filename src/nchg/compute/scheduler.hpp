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

#include <parallel_hashmap/btree.h>

#include <cstddef>
#include <filesystem>
#include <hictk/chromosome.hpp>
#include <optional>
#include <utility>

#include "nchg/expected_values.hpp"
#include "nchg/file_store.hpp"
#include "nchg/genomic_domains.hpp"
#include "nchg/tools/config.hpp"

namespace nchg {

[[nodiscard]] std::size_t process_chromosome_pair(
    const ComputePvalConfig &c, const std::optional<GenomicDomains> &domains,
    const std::optional<ExpectedValues> &expected_values = {});

[[nodiscard]] std::size_t dispatch_queries(FileStore &file_store,
                                           const ChromosomePairs &chrom_pairs,
                                           const std::optional<GenomicDomains> &domains,
                                           const std::optional<ExpectedValues> &expected_values,
                                           const ComputePvalConfig &c);

[[nodiscard]] std::filesystem::path generate_report_file_name(
    const std::filesystem::path &output_prefix);
[[nodiscard]] std::filesystem::path generate_chrom_sizes_file_name(
    const std::filesystem::path &output_prefix);
[[nodiscard]] std::filesystem::path generate_output_file_name(
    const std::filesystem::path &output_prefix, const hictk::Chromosome &chrom1,
    const hictk::Chromosome &chrom2);

}  // namespace nchg
