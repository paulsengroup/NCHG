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
        self.requires("arrow/17.0.0#81be2aa6c49800df8cc163adf4b99e9f")
        self.requires("boost/1.86.0#cd839a2082585255010f9e82eea94c7f", force=True)
        self.requires("bshoshany-thread-pool/4.1.0#be1802a8768416a6c9b1393cf0ce5e9c")
        self.requires("catch2/3.7.1#431d772165ed0bc5adaabaa44a9f53ca")
        self.requires("cli11/2.4.2#1b431bda2fb2cd3efed633899abcd8cc")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("fast_float/6.1.5#e067b96a6271d1b4c255858ca9805bdd")  # hictk
        self.requires("fmt/11.0.2#5c7438ef4d5d69ab106a41e460ce11f3", force=True)
        self.requires("hdf5/1.14.5#51799cda2ba7acaa74c9651dea284ac4", force=True)
        # self.requires("hictk/0.0.12#8e413cd45528da38b5a41ccffee41d6d")
        self.requires("highfive/2.10.0#3d1bd25944a57fa1bc30a0a22923d528")
        self.requires("libdeflate/1.22#f95aebe763153ccbc4cc76c023e42e5a")  # hictk
        self.requires("parallel-hashmap/1.37#ad1feb0b318d094645b696b7cbb0657e")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")  # hictk
        self.requires("spdlog/1.14.1#972bbf70be1da4bc57ea589af0efde03")
        self.requires("thrift/0.20.0#560fdab2e1636d4d8a0556fcf6470b89", force=True)
        self.requires("zstd/1.5.6#afefe79a309bc2a7b9f56c2093504c8b", force=True)

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
        self.options["boost"].without_charconv = True
        self.options["boost"].without_chrono = True
        self.options["boost"].without_cobalt = True
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
        self.options["zstd"].build_programs = False
