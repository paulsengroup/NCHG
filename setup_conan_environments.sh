#!/usr/bin/env bash

# Copyright (C) 2024 Roberto Rossini (roberros@uio.no)
# SPDX-License-Identifier: GPL-3.0
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public
# License along with this library.  If not, see
# <https://www.gnu.org/licenses/>.

set -e
set -u

# shellcheck disable=SC2064
trap "cd '$PWD'" EXIT

git_root="$(readlink -f "$(git rev-parse --show-toplevel)")"

wd="$git_root/conan-envs"
conanfile="$git_root/conanfile.txt"

for compiler in gcc clang; do
  for build_type in Debug Release RelWithDebInfo; do
    CC="$compiler"
    if [[ "$compiler" == gcc* ]]; then
      CXX="${compiler/gcc/g++}${compiler#gcc}"
      profile=gcc
    else
      CXX="${compiler/clang/clang++}${compiler#clang}"
      profile=clang
    fi

    export CC
    export CXX

    outdir="$wd/$compiler/$build_type"
    rm -rf "$outdir"
    mkdir -p "$outdir"

    conan install \
      --build=missing \
      --update \
      -pr "$profile"  \
      -s compiler.cppstd=17 \
      -s build_type="$build_type" \
      --output-folder="$outdir" \
      "$conanfile"

     conan install \
       --build=missing \
       --update \
       -pr "$profile"  \
       -s compiler.cppstd=17 \
       -s build_type="$build_type" \
       -o shared=True \
       --output-folder="$outdir" \
       "$conanfile"
  done
done
