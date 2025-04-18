# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: GPL-3.0
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public License as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public License along
# with this library.  If not, see
# <https://www.gnu.org/licenses/>.

from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"


class NCHGConan(ConanFile):
    name = "NCHG"
    description = "TODO"
    license = "GPL-3.0"
    topics = ("NCHG", "bioinformatics")
    homepage = "https://github.com/paulsengroup/NCHG"
    url = "https://github.com/paulsengroup/NCHG"
    package_type = "application"
    settings = "os", "arch", "compiler", "build_type"

    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }

    default_options = {
        "shared": False,
        "fPIC": True,
    }

    generators = "CMakeDeps"

    @property
    def _min_cppstd(self):
        return 23

    def requirements(self):
        self.requires("arrow/19.0.1#f6937fd566ecbec1eab37b40e292dfec")
        self.requires("boost/1.87.0#53c53f3d6eeb9db4a3d68573596db0e7", force=True)
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("catch2/3.8.0#2c87b60d2c85f3c8509bb209f37cbf67")
        self.requires("cli11/2.5.0#1b7c81ea2bff6279eb2150bbe06a200a")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("fast_float/8.0.0#edda0315516b2f1e7835972fdf5fc5ca")  # hictk
        self.requires("fmt/11.1.4#1fb24f082fabe20d28606d615ba93dfb", force=True)
        self.requires("glaze/5.0.2#1f0dcbd3b22dd5f2385ebdf3a252ec3f")
        self.requires("hdf5/1.14.5#51799cda2ba7acaa74c9651dea284ac4", force=True)
        self.requires("highfive/2.10.0#c975a16d7fe3655c173f8a9aab16b416")
        self.requires("libdeflate/1.23#4994bea7cf7e93789da161fac8e26a53")  # hictk
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")  # hictk
        self.requires("spdlog/1.15.1#92e99f07f134481bce4b70c1a41060e7")
        self.requires("thrift/0.20.0#560fdab2e1636d4d8a0556fcf6470b89", force=True)
        self.requires("xxhash/0.8.3#681d36a0a6111fc56e5e45ea182c19cc")
        self.requires("zstd/1.5.7#fde461c0d847a22f16d3066774f61b11", force=True)

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"] and self.settings.os == "Linux":
            self.settings.compiler.libcxx = "libstdc++11"

        self.options["arrow"].parquet = True
        self.options["arrow"].with_boost = True
        self.options["arrow"].with_thrift = True
        self.options["arrow"].with_lz4 = True
        self.options["arrow"].with_zstd = True
        self.options["boost"].system_no_deprecated = True
        self.options["boost"].asio_no_deprecated = True
        self.options["boost"].filesystem_no_deprecated = True
        self.options["boost"].filesystem_version = 4
        self.options["boost"].zlib = False
        self.options["boost"].bzip2 = False
        self.options["boost"].lzma = False
        self.options["boost"].zstd = False
        self.options["boost"].without_atomic = False
        self.options["boost"].without_charconv = True
        self.options["boost"].without_chrono = True
        self.options["boost"].without_cobalt = True
        self.options["boost"].without_container = True
        self.options["boost"].without_context = False
        self.options["boost"].without_contract = True
        self.options["boost"].without_coroutine = True
        # without_date_time is set to False to workaround https://github.com/conan-io/conan-center-index/issues/26890
        self.options["boost"].without_date_time = False
        self.options["boost"].without_exception = False
        self.options["boost"].without_fiber = True
        self.options["boost"].without_filesystem = False
        self.options["boost"].without_graph = True
        self.options["boost"].without_graph_parallel = True
        self.options["boost"].without_iostreams = True
        self.options["boost"].without_json = True
        self.options["boost"].without_locale = True
        self.options["boost"].without_log = True
        self.options["boost"].without_math = True
        self.options["boost"].without_mpi = True
        self.options["boost"].without_nowide = True
        self.options["boost"].without_process = False
        self.options["boost"].without_program_options = True
        self.options["boost"].without_python = True
        self.options["boost"].without_random = True
        self.options["boost"].without_regex = True
        self.options["boost"].without_serialization = True
        self.options["boost"].without_stacktrace = True
        self.options["boost"].without_system = False
        self.options["boost"].without_test = True
        self.options["boost"].without_thread = True
        self.options["boost"].without_timer = True
        self.options["boost"].without_type_erasure = True
        self.options["boost"].without_url = True
        self.options["boost"].without_wave = True
        self.options["fmt"].header_only = True
        self.options["hdf5"].enable_cxx = False
        self.options["hdf5"].hl = False
        self.options["hdf5"].threadsafe = False
        self.options["hdf5"].parallel = False
        self.options["hictk"].with_eigen = False
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False
        self.options["spdlog"].header_only = True
        self.options["thrift"].with_libevent = False
        self.options["thrift"].with_openssl = False
        self.options["xxhash"].utility = False
        self.options["zstd"].build_programs = False
