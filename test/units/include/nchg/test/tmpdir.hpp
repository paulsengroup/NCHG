// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <hictk/tmpdir.hpp>

namespace nchg::test {

// NOLINTBEGIN(cert-err58-cpp)
inline const hictk::internal::TmpDir testdir{true};
inline const std::filesystem::path datadir{"test/data"};
// NOLINTEND(cert-err58-cpp)

}  // namespace nchg::test
