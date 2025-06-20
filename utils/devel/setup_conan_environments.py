#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini (roberros@uio.no)
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

import argparse
import itertools
import os
import pathlib
import shlex
import shutil
import subprocess as sp
from typing import Dict


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        help="Path where to store the resulting *.cmake files",
    )
    cli.add_argument(
        "--profile",
        nargs="+",
        type=str,
        default=("gcc", "clang"),
        choices={"gcc", "clang"},
        help="Names of the conan profiles to be used.",
    )
    cli.add_argument(
        "--build-type",
        nargs="+",
        type=str,
        default=("Debug", "RelWithDebInfo", "Release"),
        help="Conan build types.",
    )
    cli.add_argument(
        "--build-shared-only",
        action="store_true",
        default=False,
        help="Build dependencies as shared libraries only.",
    )
    cli.add_argument(
        "--build-static-only",
        action="store_true",
        default=False,
        help="Build dependencies as static libraries only.",
    )
    cli.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print the commands that would be executed if --dry-run was not specified, then exit.",
    )

    return cli


def infer_root_dir() -> pathlib.Path:
    path = pathlib.Path(sp.check_output(["git", "rev-parse", "--show-toplevel"], encoding="utf-8").strip())

    if not path.is_dir():
        raise RuntimeError("Unable to infer repository root!")

    return path


def run_or_print(args, env: Dict[str, str], dry_run: bool):
    if dry_run:
        print(shlex.join(str(x) for x in args))
    else:
        sp.check_call(args, env=env)


def run_conan(
    profile: str,
    build_type: str,
    shared: bool,
    output_prefix: pathlib.Path,
    dry_run: bool,
):
    output_folder = output_prefix / profile / build_type / ("shared" if shared else "static")

    base_args = [
        "conan",
        "install",
        infer_root_dir() / "conanfile.py",
        "--build=missing",
        "--update",
        "--profile",
        profile,
        "--settings",
        f"build_type={build_type}",
        "--options",
        f"*/*:shared={shared}",
        "--output-folder",
        output_folder,
    ]

    env = os.environ.copy()
    env["CC"] = profile
    if profile == "gcc":
        env["CXX"] = "g++"
    elif profile == "clang":
        env["CXX"] = "clang++"
    else:
        raise RuntimeError(f'Unrecognized compiler "{profile}". Profiles should be either named "gcc" or "clang"')

    run_or_print(base_args + ["--settings", "compiler.cppstd=23"], env=env, dry_run=dry_run)
    run_or_print(
        base_args + ["--settings", "compiler.cppstd=23", "--options", "NCHG/*:with_glaze_only=True"],
        env=env,
        dry_run=dry_run,
    )
    run_or_print(
        base_args + ["--settings", "compiler.cppstd=20", "--options", "NCHG/*:with_duckdb_only=True"],
        env=env,
        dry_run=dry_run,
    )


def main():
    args = vars(make_cli().parse_args())

    profiles = args["profile"]
    build_types = args["build_type"]
    if args["build_shared_only"]:
        shared_build = [True]
    elif args["build_static_only"]:
        shared_build = [False]
    else:
        shared_build = [True, False]

    dry_run = args["dry_run"]

    output_prefix = infer_root_dir() / "conan-envs"
    if output_prefix.exists():
        shutil.rmtree(output_prefix)

    output_prefix.mkdir(exist_ok=True)

    for args in itertools.product(profiles, build_types, shared_build):
        run_conan(*args, output_prefix=output_prefix, dry_run=dry_run)


if __name__ == "__main__":
    main()
