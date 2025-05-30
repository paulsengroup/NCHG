#!/usr/bin/env bash

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

set -eu
set -o pipefail

if [ $# -ne 1 ]; then
  2>&1 echo "Usage: $0 nchg:latest"
  exit 1
fi

IMG="$1"

tmpdir="$(mktemp -d)"
# shellcheck disable=SC2064
trap "rm -rf '$tmpdir'" EXIT

cat > "$tmpdir/runme.sh" <<- 'EOM'

set -eu

whereis -b NCHG
NCHG --help
NCHG --version

if NCHG compute \
  /tmp/data/ENCFF447ERX.minified.mcool \
  --resolution 1000000 \
  --threads "$(nproc)" \
  out/ENCFF447ERX.1000000; then
  2>&1 echo "### SUCCESS ###"
else
  2>&1 echo "### FAILURE ###"
  exit 1
fi

EOM

chmod 755 "$tmpdir/runme.sh"

sudo docker run --rm --entrypoint=/bin/bash \
  -v "$tmpdir/runme.sh:/tmp/runme.sh:ro" \
  -v "$PWD/test/scripts:/tmp/nghc/test/scripts:ro" \
  -v "$PWD/test/data/:/tmp/data/:ro" \
  "$IMG" \
  /tmp/runme.sh
