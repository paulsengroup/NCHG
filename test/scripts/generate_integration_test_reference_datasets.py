#!/usr/bin/env python3

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


import argparse
import logging
import math
import pathlib
import shutil
import subprocess as sp
import sys
import tempfile
from typing import Any, Dict, List

import pandas as pd


def existing_file(path: str) -> pathlib.Path:
    p = pathlib.Path(path)
    if p.is_file():
        return p

    raise RuntimeError(f'"{path}" does not point to an existing file')


def duration(s: str) -> float:
    try:
        n = float(s)
        if not math.isnan(n) and n > 0:
            return n
    except Exception:  # noqa
        pass

    raise RuntimeError(f'"{s}" is not a valid duration')


def non_zero_int(s: str) -> int:
    try:
        n = int(s)
        if n > 0:
            return n
    except Exception:  # noqa
        pass

    raise RuntimeError(f'"{s}" is not a valid non-zero integer')


def executable_file(path: str) -> pathlib.Path:
    exe = shutil.which(path)
    if exe:
        return pathlib.Path(exe)

    raise RuntimeError(f'Unable to find executable "{path}"')


def add_common_flags(parser):
    parser.add_argument(
        "--nchg-bin",
        type=executable_file,
        help="Path to NCHG's executable.",
    )

    parser.add_argument(
        "--nproc",
        type=non_zero_int,
        default=1,
        help="Maximum number of parallel processes to spawn.",
    )

    parser.add_argument(
        "--verbosity",
        choices={"debug", "info", "warnings", "error", "critical"},
        type=str,
        default="info",
        help="Tweak log verbosity.",
    )

    parser.add_argument(
        "--timeout",
        type=duration,
        default=60,
        help="Subprocess timeout (seconds).",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing file(s).",
    )


def make_cross_product_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "cross-product",
        help="Generate the reference dataset for NCHG cross-product.",
    )

    sc.add_argument(
        "1d-domains",
        type=existing_file,
        help="Path to a BED3+ files with the list of domains to be processed.",
    )

    sc.add_argument(
        "out-path",
        type=pathlib.Path,
        help="Output path where to store the resulting domains.",
    )

    sc.add_argument(
        "--type",
        type=str,
        choices={"gw", "cis", "trans"},
        default="gw",
        help="Type of domains to be generated.",
    )

    sc.add_argument(
        "--verbosity",
        choices={"debug", "info", "warnings", "error", "critical"},
        type=str,
        default="info",
        help="Tweak log verbosity.",
    )

    sc.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing file(s).",
    )


def make_compute_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "compute",
        help="Generate the reference dataset for NCHG compute.",
    )

    sc.add_argument(
        "cool-uri",
        type=str,
        help="Path or URI to an existing .cool file.",
    )
    sc.add_argument(
        "domains",
        type=existing_file,
        help="Path to a BED3+ file with the list of domains to be processed.",
    )
    sc.add_argument(
        "out-prefix",
        type=pathlib.Path,
        help="Output path to pass to NCHG compute.",
    )

    add_common_flags(sc)


def make_merge_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "merge",
        help="Generate the reference dataset for NCHG merge.",
    )

    sc.add_argument(
        "input-prefix",
        type=pathlib.Path,
        help="Path prefix to pass to NCHG merge.",
    )
    sc.add_argument(
        "output-parquet",
        type=pathlib.Path,
        help="Path where to store NCHG merge's output.",
    )

    add_common_flags(sc)


def make_filter_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "filter",
        help="Generate the reference dataset for NCHG filter.",
    )

    sc.add_argument(
        "input-parquet",
        type=existing_file,
        help="Path to the .parquet file to be passed to NCHG filter.",
    )
    sc.add_argument(
        "output-parquet",
        type=pathlib.Path,
        help="Path where to store NCHG filter's output.",
    )

    add_common_flags(sc)


def make_view_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "view",
        help="Generate the reference dataset for NCHG view.",
    )

    sc.add_argument(
        "input-parquet",
        type=existing_file,
        help="Path to the .parquet file to be passed to NCHG view.",
    )
    sc.add_argument(
        "output-tsv",
        type=pathlib.Path,
        help="Path where to store NCHG view's output.",
    )

    add_common_flags(sc)


def make_expected_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "expected",
        help="Generate the reference dataset for NCHG expected.",
    )

    sc.add_argument(
        "cool-uri",
        type=str,
        help="Path or URI to an existing .cool file.",
    )
    sc.add_argument(
        "output-h5",
        type=pathlib.Path,
        help="Path where to store NCHG expect's output.",
    )
    sc.add_argument(
        "--masked-intervals",
        nargs="+",
        type=str,
        help="One or more intervals to be masked out (UCSC format).",
    )

    add_common_flags(sc)


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser("Generate the reference datasets for NCHG integration tests.")

    cli.add_argument(
        "cool-uri",
        type=str,
        help="Path or URI to an existing .cool file.",
    )
    cli.add_argument(
        "output-dir",
        type=pathlib.Path,
        help="Output folder where to store the output produced by NCHG.",
    )
    cli.add_argument(
        "--domains",
        type=existing_file,
        help="Path to a BED3+ file with the list of domains to be processed.",
    )
    add_common_flags(cli)

    sub_parser = cli.add_subparsers(
        title="subcommands",
        dest="command",
        help="List of available subcommands:",
    )

    make_compute_sc(sub_parser)
    make_merge_sc(sub_parser)
    make_filter_sc(sub_parser)
    make_view_sc(sub_parser)
    make_expected_sc(sub_parser)

    return cli


def run_nchg_command(nchg_bin: pathlib.Path, subcmd: str, *args, **kwargs):
    if "timeout" not in kwargs:
        kwargs["timeout"] = 60
    else:
        assert kwargs["timeout"] > 0

    if "force" not in kwargs:
        kwargs["force"] = False

    cmd = [nchg_bin, subcmd]
    cmd.extend((str(x) for x in args))

    if kwargs["force"]:
        cmd.append("--force")

    logging.debug("launching %s...", [str(x) for x in cmd])

    res = sp.run(
        cmd,
        stdin=sp.DEVNULL,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        timeout=kwargs["timeout"],
        encoding="utf-8",
    )

    for line in res.stdout.strip().split("\n"):
        logging.info("NCHG %s [stdout]: %s", subcmd, line)
    for line in res.stderr.strip().split("\n"):
        logging.info("NCHG %s [stderr]: %s", subcmd, line)

    if res.returncode != 0:
        raise RuntimeError(f"{cmd} returned with exit code {res.returncode}")


def generate_2d_domains(
    domains: pathlib.Path,
    output_path: pathlib.Path,
    domain_type: str,
    force: bool,
) -> pathlib.Path:
    if output_path.exists() and not force:
        raise RuntimeError(f'refusing to overwrite file "{output_path}": pass --force to overwrite.')

    parent_dir = output_path.parent
    parent_dir.mkdir(exist_ok=True)

    logging.info("reading domains from %s...", domains)

    df = pd.read_table(domains, names=["chrom", "start", "end"], usecols=[0, 1, 2])
    chroms = {name: i for i, name in enumerate(df["chrom"].unique())}

    logging.info("parsed %d chromosomes...", len(chroms))
    df["id"] = df["chrom"].map(chroms)

    df = df.merge(df, how="cross", suffixes=("1", "2"))
    # Drop domain pairs overlapping with the lower-triangular matrix
    df = df[(df["id1"] < df["id2"]) | ((df["id1"] == df["id2"]) & (df["start1"] <= df["start2"]))]

    if domain_type != "cis":
        df = df.sort_values(["id1", "id2", "start1", "start2", "end1", "end2"])

    df = df.reset_index(drop=True).drop(columns=["id1", "id2"])

    size = len(df)
    logging.info("generated %d domain pairs...", size)
    if domain_type == "cis":
        logging.info("dropping trans domains...")
        df = df[df["chrom1"] == df["chrom2"]]
        logging.info("dropped %d domains...", size - len(df))
    elif domain_type == "trans":
        logging.info("dropping cis domains...")
        df = df[df["chrom1"] != df["chrom2"]]
        logging.info("dropped %d domains...", size - len(df))

    logging.info("writing 2D domains to %s...", output_path)
    df.to_csv(output_path, sep="\t", index=False, header=False)

    return output_path


def run_nchg_compute(
    nchg_bin: pathlib.Path,
    uri: str,
    out_prefix: pathlib.Path,
    domains: pathlib.Path | None,
    nproc: int,
    force: bool,
    timeout: float,
) -> pathlib.Path:
    args = [uri, out_prefix, "--threads", nproc]
    if domains is not None:
        args.extend(("--domains", domains))

    run_nchg_command(
        nchg_bin,
        "compute",
        *args,
        force=force,
        timeout=timeout,
    )
    return out_prefix


def run_nchg_merge(
    nchg_bin: pathlib.Path, input_prefix: pathlib.Path, output_path: pathlib.Path, force: bool, timeout: float
) -> pathlib.Path:
    run_nchg_command(
        nchg_bin,
        "merge",
        input_prefix,
        output_path,
        force=force,
        timeout=timeout,
    )
    assert output_path.is_file()
    return output_path


def run_nchg_filter(
    nchg_bin: pathlib.Path, input_parquet: pathlib.Path, output_parquet: pathlib.Path, force: bool, timeout: float
) -> pathlib.Path:
    run_nchg_command(
        nchg_bin,
        "filter",
        input_parquet,
        output_parquet,
        force=force,
        timeout=timeout,
    )
    assert output_parquet.is_file()
    return output_parquet


def run_nchg_view(
    nchg_bin: pathlib.Path, input_parquet: pathlib.Path, output_tsv: pathlib.Path, force: bool, timeout: float
) -> pathlib.Path:
    assert timeout > 0

    cmd = [nchg_bin, "view", input_parquet]
    logging.debug("launching %s...", [str(x) for x in cmd])

    if output_tsv.exists() and not force:
        raise RuntimeError(f'refusing to overwrite file "{output_tsv}": pass --force to overwrite.')

    with output_tsv.open("w") as f:
        res = sp.run(
            cmd,
            stdin=sp.DEVNULL,
            stdout=f,
            stderr=sp.PIPE,
            timeout=timeout,
        )
        for line in res.stderr:
            logging.info("NCHG view [stderr]: %s", line)

        if res.returncode != 0:
            output_tsv.unlink(missing_ok=True)
            raise RuntimeError(f"{cmd} returned with exit code {res.returncode}")

    return output_tsv


def run_nchg_expected(
    nchg_bin: pathlib.Path,
    uri: str,
    output_h5: pathlib.Path,
    masked_regions: List[str] | None,
    force: bool,
    timeout: float,
) -> pathlib.Path:
    with tempfile.NamedTemporaryFile() as tmpfile:
        for region in masked_regions:
            tmpfile.write(f"{region}\n")

        args = [nchg_bin, "expected", uri, "--output", output_h5]
        if len(masked_regions) != 0:
            tmpfile.flush()
            args.extend(("--bin-mask", tmpfile.name))

        run_nchg_command(
            *args,
            force=force,
            timeout=timeout,
        )

    assert output_h5.exists()
    return output_h5


def generate_all_files(nchg_bin: pathlib.Path, force: bool, timeout: float, args: Dict[str, Any]):
    suffix = pathlib.Path(args["cool-uri"].partition("::")[0]).stem

    for dom_type in ["gw", "cis", "trans"]:
        output = args["output-dir"] / "cross_product" / f"{suffix}.{dom_type}-domains.bedpe"
        generate_2d_domains(args["domains"], output, dom_type, args["force"])

    domains = args["output-dir"] / "cross_product" / f"{suffix}.gw-domains.bedpe"
    outprefix = args["output-dir"] / "compute_with_domains" / suffix
    run_nchg_compute(
        nchg_bin,
        args["cool-uri"],
        outprefix,
        domains,
        args["nproc"],
        force,
        timeout,
    )

    outprefix = args["output-dir"] / "compute" / suffix
    run_nchg_compute(
        nchg_bin,
        args["cool-uri"],
        outprefix,
        None,
        args["nproc"],
        force,
        timeout,
    )

    output = args["output-dir"] / "merge" / f"{suffix}.parquet"
    output.parent.mkdir(parents=True, exist_ok=True)
    run_nchg_merge(
        nchg_bin,
        outprefix,
        output,
        force,
        timeout,
    )

    input = output
    output = args["output-dir"] / "filter" / f"{suffix}.parquet"
    output.parent.mkdir(parents=True, exist_ok=True)
    run_nchg_filter(
        nchg_bin,
        input,
        output,
        force,
        timeout,
    )

    input = output
    output = args["output-dir"] / "view" / f"{suffix}.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    run_nchg_view(
        nchg_bin,
        input,
        output,
        force,
        timeout,
    )

    output = args["output-dir"] / "expected" / f"{suffix}.h5"
    output.parent.mkdir(parents=True, exist_ok=True)
    run_nchg_expected(
        nchg_bin,
        args["cool-uri"],
        output,
        [],
        force,
        timeout,
    )


def find_nchg_exec(path: pathlib.Path | None) -> pathlib.Path:
    if path:
        if shutil.which(path):
            return path
        raise RuntimeError(f'"{path}" does not seem to be a valid executable.')

    path = shutil.which("NCHG")
    if path is not None:
        return pathlib.Path(path)

    raise RuntimeError(
        "unable to find NCHG in your PATH: please specify the path to NCHG's executable with --nchg-bin."
    )


def main():
    args = vars(make_cli().parse_args())

    setup_logger(args["verbosity"].upper())

    nchg_bin = find_nchg_exec(args["nchg_bin"])
    force = args["force"]
    timeout = args["timeout"]

    cmd = args["command"]

    if cmd == "cross-product":
        generate_2d_domains(
            args["domains"],
            args["out-path"],
            args["type"],
            args["force"],
        )
    elif cmd == "compute":
        run_nchg_compute(
            nchg_bin,
            args["cool-uri"],
            args["out-prefix"],
            args["domains"],
            args["nproc"],
            force,
            timeout,
        )
    elif cmd == "merge":
        run_nchg_merge(
            nchg_bin,
            args["input-prefix"],
            args["output-parquet"],
            force,
            timeout,
        )
    elif cmd == "filter":
        run_nchg_filter(
            nchg_bin,
            args["input-parquet"],
            args["output-parquet"],
            force,
            timeout,
        )
    elif cmd == "view":
        run_nchg_view(
            nchg_bin,
            args["input-parquet"],
            args["output-tsv"],
            force,
            timeout,
        )
    elif cmd == "expected":
        run_nchg_expected(
            nchg_bin,
            args["cool-uri"],
            args["output-h5"],
            args["masked_intervals"],
            force,
            timeout,
        )
    elif not cmd:
        generate_all_files(
            nchg_bin,
            force,
            timeout,
            args,
        )
    else:
        raise NotImplementedError


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    sys.exit(main())
