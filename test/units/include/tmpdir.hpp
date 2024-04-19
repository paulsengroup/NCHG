// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <hictk/tmpdir.hpp>

namespace nchg::test {

inline const hictk::internal::TmpDir testdir{true};       // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)

}  // namespace nchg::test
