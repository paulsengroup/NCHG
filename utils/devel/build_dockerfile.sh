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

set -e
set -u
set -o pipefail

IMAGE_NAME='nchg'

GIT_HASH="$(git rev-parse HEAD)"
GIT_SHORT_HASH="$(git rev-parse --short HEAD)"
GIT_TAG="$(git for-each-ref 'refs/tags/v*.*.*' --count 1 --sort=-v:refname --format "%(refname:short)"  --points-at HEAD)"
CREATION_DATE="$(date --iso-8601)"

if [[ $(git status --porcelain -uno) ]]; then
  GIT_IS_DIRTY=1
else
  GIT_IS_DIRTY=0
fi

IMAGE_TAG="sha-$GIT_SHORT_HASH"
if [ $GIT_IS_DIRTY -ne 0 ]; then
  IMAGE_TAG+='-dirty'
fi

if [ -z "$GIT_TAG" ]; then
  GIT_TAG="sha-$GIT_SHORT_HASH"
else
  GIT_TAG="$GIT_TAG"
fi

2>&1 echo "Building \"$IMAGE_NAME:$IMAGE_TAG\"..."

sudo docker pull docker.io/library/ubuntu:22.04
FINAL_BASE_IMAGE_DIGEST="$(sudo docker inspect --format='{{index .RepoDigests 0}}' docker.io/library/ubuntu:22.04 | grep -o '[[:alnum:]:]\+$')"

BUILD_BASE_IMAGE='ghcr.io/paulsengroup/ci-docker-images/ubuntu-22.04-cxx-clang-18:latest'

sudo docker pull "$BUILD_BASE_IMAGE"

# sudo docker buildx build --platform linux/amd64,linux/arm64 \
sudo docker buildx build --platform linux/amd64 \
  --build-arg "BUILD_BASE_IMAGE=$BUILD_BASE_IMAGE" \
  --build-arg "FINAL_BASE_IMAGE=docker.io/library/ubuntu" \
  --build-arg "FINAL_BASE_IMAGE_TAG=22.04" \
  --build-arg "FINAL_BASE_IMAGE_DIGEST=$FINAL_BASE_IMAGE_DIGEST" \
  --build-arg "C_COMPILER=clang-18" \
  --build-arg "CXX_COMPILER=clang++-18" \
  --build-arg "GIT_HASH=$GIT_HASH" \
  --build-arg "GIT_SHORT_HASH=$GIT_SHORT_HASH" \
  --build-arg "GIT_TAG=$GIT_TAG" \
  --build-arg "GIT_IS_DIRTY=$GIT_IS_DIRTY" \
  --build-arg "CREATION_DATE=$CREATION_DATE" \
  -t "$IMAGE_NAME:latest" \
  -t "$IMAGE_NAME:$(echo "$CREATION_DATE" | tr -d '\-' )" \
  -t "$IMAGE_NAME:$IMAGE_TAG" \
  "$(git rev-parse --show-toplevel)"

 # sudo singularity build -F "${img_name}_v${ver}.sif" \
 #                           "docker-daemon://${img_name}:${ver}"
