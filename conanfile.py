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
        return 17

    def requirements(self):
        self.requires("arrow/16.0.0#c0b940c80fddf319983165f7588c279d")
        self.requires("boost/1.84.0#5cc9767cd47ba6b5e1b53ab4670f07fc", override=True)
        self.requires("bshoshany-thread-pool/4.1.0#be1802a8768416a6c9b1393cf0ce5e9c")
        self.requires("catch2/3.5.4#d346ca291f8f62040fd9c1a891654711")
        self.requires("cli11/2.4.1#afacffd31f631bbb8b7c7d6425fe7a66")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("fmt/10.2.1#9199a7a0611866dea5c8849a77467b25")
        self.requires("hdf5/1.14.3#31ccd8d4de83844f5db48471df1944a1")
        self.requires("hictk/0.0.12#8e413cd45528da38b5a41ccffee41d6d")
        self.requires("highfive/2.9.0#c57477beed8b0110fadeb6da8f48bcc5")
        self.requires("parallel-hashmap/1.3.11#1e67f4855a3f7cdeb977cc472113baf7")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("spdlog/1.13.0#8e88198fd5b9ee31d329431a6d0ccaa2")
        self.requires("thrift/0.18.1#4e5674c24f99dde562c3926f9cb2ff9d", override=True)
        self.requires("zstd/1.5.6#67383dae85d33f43823e7751a6745ea1", override=True)

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"]:
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
        self.options["boost"].without_chrono = True
        self.options["boost"].without_container = True
        self.options["boost"].without_context = True
        self.options["boost"].without_contract = True
        self.options["boost"].without_coroutine = True
        self.options["boost"].without_date_time = True
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
        self.options["zstd"].build_programs = False
