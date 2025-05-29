#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT


import argparse
import logging
import multiprocessing as mp
import pathlib
from typing import List

import hictkpy


def existing_file(arg: str) -> pathlib.Path:
    if (path := pathlib.Path(arg)).is_file():
        return path

    raise argparse.ArgumentTypeError(f'Not an existing file: "{arg}"')


def positive_int(arg) -> int:
    if (n := int(arg)) > 0:
        return n

    raise argparse.ArgumentTypeError("Not a positive int")


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "hic-matrix",
        type=existing_file,
        help="Path to the Hi-C matrix in .hic or .[m]cool format to be processed.",
    )
    cli.add_argument(
        "--chromosomes",
        required=True,
        type=str,
        nargs="+",
        help="Name of one or more chromosomes to be processed. "
        "The resulting matrix will use the same reference genome as the input matrix, but only contain interactions for the given chromosome(s).",
    )
    cli.add_argument(
        "--resolution",
        type=positive_int,
        help="Resolution of the Hi-C matrix in bp. Required when input matrix is multi-resolution.",
    )
    cli.add_argument(
        "--output-file",
        required=True,
        type=pathlib.Path,
        help="Path to the output file to be created. "
        "File format is inferred from the extension, which should be either .cool or .hic.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing file(s).",
    )

    return cli


def import_chromosomes(chroms: List[str], f: hictkpy.File) -> List[str]:
    indexed_chroms = {chrom: i for i, chrom in enumerate(f.chromosomes())}
    valid_chroms = {}

    for chrom in chroms:
        if chrom in indexed_chroms:
            valid_chroms[chrom] = indexed_chroms[chrom]
        else:
            raise RuntimeError(f'Chromosome "{chrom}" not found in file {f.uri()}')

    # Return chromosomes sorted in by the order they appear in the reference genome
    return [chrom for chrom, _ in sorted(valid_chroms.items(), key=lambda x: x[1])]


def create_output_file(
    path: pathlib.Path,
    input_matrix: hictkpy.File,
    force: bool,
) -> hictkpy.cooler.FileWriter | hictkpy.hic.FileWriter:
    if not force and path.exists():
        raise RuntimeError(f'Refusing to overwrite file "{path}". Pass --force to overwrite.')

    path.unlink(missing_ok=True)

    if input_matrix.bins().type() != "fixed":
        raise RuntimeError(f"Files with bin-type={input_matrix.bins().type()} are not supported.")

    if path.suffix == ".cool":
        logging.info('initializing output file "%s" using Cooler format', path)
        return hictkpy.cooler.FileWriter(
            path,
            input_matrix.chromosomes(),
            input_matrix.resolution(),
            input_matrix.attributes().get("assembly", "unknown"),
            compression_lvl=9,
        )

    if path.suffix == ".hic":
        logging.info('initializing output file "%s" using .hic format', path)
        return hictkpy.hic.FileWriter(
            path,
            input_matrix.chromosomes(),
            input_matrix.resolution(),
            input_matrix.attributes().get("assembly", "unknown"),
            compression_lvl=12,
            n_threads=mp.cpu_count(),
        )

    raise RuntimeError(f"Unable to infer file format from path with extension {path.suffix}.")


def setup_logger(level: str = "INFO"):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    f = hictkpy.File(args["hic-matrix"], args["resolution"])
    chroms = import_chromosomes(args["chromosomes"], f)

    logging.info("will process interactions for the following %d chromosome(s): %s", len(chroms), ", ".join(chroms))

    writer = create_output_file(pathlib.Path(args["output_file"]), f, args["force"])

    for i, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i:]:
            logging.info("processing interactions for %s:%s matrix", chrom1, chrom2)
            df = f.fetch(chrom1, chrom2).to_pandas()
            writer.add_pixels(df)

    writer.finalize(log_lvl="info")

    dest = writer.path()
    if hictkpy.is_cooler(dest):
        dest = writer.path().with_suffix(".mcool")

    print("\n")
    print("# To finish generating the test dataset, please run the following commands:")
    print(f"mv {writer.path()} {writer.path()}.tmp")
    print(f"hictk zoomify {writer.path()}.tmp {dest} -l 12")
    print(f"hictk balance vc {dest} --mode=cis --no-create-weight-link")
    print(f"hictk balance scale {dest} --mode=cis --threads={mp.cpu_count()} --no-create-weight-link  # --in-memory")
    print(f"hictk balance ice {dest} --mode=cis --threads={mp.cpu_count()} --create-weight-link  # --in-memory")
    print(f"# rm {writer.path()}.tmp  # optional")


if __name__ == "__main__":
    setup_logger()
    main()
