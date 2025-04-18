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
import gzip
import logging
import pathlib
import shutil
import subprocess as sp
import sys
import tempfile
from urllib.request import urlretrieve

import bioframe as bf
import cooler
import cooltools
import pandas as pd


def existing_file(path: str) -> pathlib.Path:
    p = pathlib.Path(path)
    if p.is_file():
        return p

    raise RuntimeError(f'"{path}" does not point to an existing file')


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "cooler",
        type=str,
        help="Path to a cooler file (URI syntax supported).",
    )
    cli.add_argument(
        "--ref-genome-fna",
        type=existing_file,
        help="Path to a FASTA file with the reference genome for compartment phasing.\n"
        "When not provided, the hg38 reference genome is downloaded from the UCSC FTP server.",
    )
    cli.add_argument(
        "--verbosity",
        choices={"debug", "info", "warnings", "error", "critical"},
        type=str,
        default="info",
        help="Tweak log verbosity.",
    )

    return cli


def read_ref_genome(path: pathlib.Path) -> pd.DataFrame:
    logging.info("reading reference genome from %s...", path)
    return bf.load_fasta(path.as_posix())


def decompress_gz(path: pathlib.Path):
    logging.info("decompressing file %s...", path)
    dest = path.with_suffix(".tmp")
    if shutil.which("gzip"):
        with dest.open("wb") as f:
            sp.run(["gzip", "-dc", path], stdin=sp.DEVNULL, stdout=f)
    else:
        with gzip.open(path, "rb") as f1, dest.open("wb") as f2:
            shutil.copyfileobj(f1, f2)

    path.unlink()
    return dest.rename(path)


def download_hg38_fna(tmpdir: pathlib.Path) -> pd.DataFrame:
    url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"

    logging.info("downloading hg38.fa.gz from %s...", url)

    with tempfile.NamedTemporaryFile(dir=tmpdir) as tmpfile:
        path = pathlib.Path(tmpfile.name)
        urlretrieve(url, path)
        return read_ref_genome(decompress_gz(path))


def get_ref_genome(path: pathlib.Path | None, tmpdir: pathlib.Path) -> pd.DataFrame:
    if path is None:
        return download_hg38_fna(tmpdir)

    return read_ref_genome(path)


def compute_phasing_track(path: pathlib.Path | None, bins: pd.DataFrame) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        ref_genome = get_ref_genome(path, tmpdir)
        logging.info("computing phasing track...")
        return bf.frac_gc(bins, ref_genome)


def compute_eigenvect(clr: cooler.Cooler, phasing_track: pd.DataFrame) -> pd.DataFrame:
    logging.info("computing eigenvector for interactions from Cooler at %s...", clr.uri)
    cis_eigs = cooltools.eigs_cis(
        clr,
        phasing_track,
        n_eigs=1,
    )

    return cis_eigs[1][["chrom", "start", "end", "E1"]]


def cluster_domains(df: pd.DataFrame, chroms: pd.Series) -> pd.DataFrame:
    df["compartment"] = "A"
    df.loc[df["E1"] < 0, "compartment"] = "B"

    df = bf.merge(df, min_dist=0, on=["compartment"]).drop(columns=["n_intervals"])

    view = chroms.to_frame().reset_index().rename(columns={"name": "chrom", "length": "end"})
    view["start"] = 0

    return bf.sort_bedframe(df, view)


def main():
    args = vars(make_cli().parse_args())
    setup_logger(args["verbosity"].upper())

    clr = cooler.Cooler(args["cooler"])
    phasing_track = compute_phasing_track(args["ref_genome_fna"], clr.bins()[:][["chrom", "start", "end"]])

    df = compute_eigenvect(clr, phasing_track)
    df = cluster_domains(df, clr.chromsizes)

    df.to_csv(sys.stdout, sep="\t", header=False, index=False)


def setup_logger(level):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    main()
