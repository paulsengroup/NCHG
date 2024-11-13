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
import pathlib
import sys
from typing import List, Set, Tuple

import h5py
import numpy as np
import pandas as pd


def existing_file(path: str) -> pathlib.Path:
    p = pathlib.Path(path)
    if p.is_file():
        return p

    raise RuntimeError(f'"{path}" does not point to an existing file')


def add_common_flags(parser):
    parser.add_argument(
        "--verbosity",
        choices={"debug", "info", "warnings", "error", "critical"},
        type=str,
        default="info",
        help="Tweak log verbosity.",
    )


def make_cartesian_product_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "cartesian-product",
        help="Validate the output produced by NCHG cartesian-product.",
    )
    sc.add_argument(
        "test-bedpe",
        type=pathlib.Path,
        help="Path to the BEDPE file to be tested.",
    )
    sc.add_argument(
        "ref-bedpe",
        type=pathlib.Path,
        help="Path to the BEDPE file to be used as ground truth.",
    )

    add_common_flags(sc)


def make_compute_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "compute",
        help="Validate the output produced by NCHG compute.",
    )
    sc.add_argument(
        "test-prefix",
        type=pathlib.Path,
        help="Path to the prefix corresponding to the files to be tested.",
    )
    sc.add_argument(
        "ref-prefix",
        type=pathlib.Path,
        help="Path to the prefix corresponding to the files to be used as ground truth.",
    )

    add_common_flags(sc)


def make_merge_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "merge",
        help="Validate the output produced by NCHG merge.",
    )
    sc.add_argument(
        "test-parquet",
        type=existing_file,
        help="Path to the .parquet file to be tested.",
    )
    sc.add_argument(
        "ref-parquet",
        type=existing_file,
        help="Path to the .parquet file to be used as ground truth.",
    )

    add_common_flags(sc)


def make_filter_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "filter",
        help="Validate the output produced by NCHG filter.",
    )
    sc.add_argument(
        "test-parquet",
        type=existing_file,
        help="Path to the .parquet file to be tested.",
    )
    sc.add_argument(
        "ref-parquet",
        type=existing_file,
        help="Path to the .parquet file to be used as ground truth.",
    )

    add_common_flags(sc)


def make_view_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "view",
        help="Validate the output produced by NCHG view.",
    )
    sc.add_argument(
        "test-tsv",
        type=existing_file,
        help="Path to the .tsv file to be tested.",
    )
    sc.add_argument(
        "ref-tsv",
        type=existing_file,
        help="Path to the .tsv file to be used as ground truth.",
    )

    add_common_flags(sc)


def make_expected_sc(main_parser):
    sc: argparse.ArgumentParser = main_parser.add_parser(
        "expected",
        help="Validate the output produced by NCHG expected.",
    )
    sc.add_argument(
        "test-h5",
        type=existing_file,
        help="Path to the .h5 file to be tested.",
    )
    sc.add_argument(
        "ref-h5",
        type=existing_file,
        help="Path to the .h5 file to be used as ground truth.",
    )

    add_common_flags(sc)


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser(
        "Validate the output file(s) produced by NCHG.",
    )
    sub_parser = cli.add_subparsers(
        title="subcommands", dest="command", required=True, help="List of available subcommands:"
    )

    make_cartesian_product_sc(sub_parser)
    make_compute_sc(sub_parser)
    make_merge_sc(sub_parser)
    make_filter_sc(sub_parser)
    make_view_sc(sub_parser)
    make_expected_sc(sub_parser)

    return cli


def import_chrom_sizes(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_table(path, names=["chrom", "size"])


def import_table(path: pathlib.Path) -> pd.DataFrame:
    try:
        return pd.read_parquet(path)
    except Exception:  # noqa
        pass

    try:
        return pd.read_table(path)
    except Exception as e:
        raise RuntimeError(f"failed to read data from {path}: file is not in .tsv or .parquet format") from e


def import_bedpe(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_table(path, names=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])


def validate_columns(expected: pd.DataFrame, found: pd.DataFrame, path: pathlib.Path):
    logging.debug("validating columns for %s...", path)

    if len(expected.columns) != len(found.columns):
        raise RuntimeError(f"column number mismatch: expected {len(expected.columns)}, found {len(found.columns)}")

    if (expected.columns != found.columns).any():
        raise RuntimeError(
            f"column name mismatch: expected {expected.columns.tolist()}, found {found.columns.tolist()}"
        )

    if (expected.dtypes != found.dtypes).any():
        raise RuntimeError(f"column dtype mismatch: expected {expected.dtypes}, found {found.dtypes}")

    logging.debug("column validation for %s was successful", path)


def validate_chrom_sizes(expected: pd.DataFrame, found: pd.DataFrame, path: pathlib.Path):
    logging.debug("validating %s...", path)
    if len(expected) != len(found):
        raise RuntimeError(
            f"unexpected number of chromosomes: expected {len(expected)} chromosomes, found {len(expected)}."
        )

    validate_columns(expected, found, path)

    if (expected["chrom"] != found["chrom"]).any():
        raise RuntimeError(
            f"chromosome name validation failed: expected {expected["chrom"].tolist()}, found {found["chrom"].tolist()}"
        )

    if (expected["size"] != found["size"]).any():
        raise RuntimeError(
            f"chromosome size validation failed: expected {expected["size"].tolist()}, found {found["size"].tolist()}"
        )

    logging.debug("%s validation was successful!", path)


def validate_pvalue_tables(expected: pd.DataFrame, found: pd.DataFrame, path: pathlib.Path):
    logging.debug("validating %s...", path)

    if len(expected) != len(found):
        raise RuntimeError(f"table length is not correct: expected {len(expected)} rows, found {len(found)}")

    df = expected.merge(
        found,
        how="outer",
        on=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
        suffixes=("_expected", "_found"),
    )

    errors = []
    for col in ["pvalue", "observed_count", "expected_count", "log_ratio", "odds_ratio", "omega"]:
        num_mismatches = (~np.isclose(df[f"{col}_expected"], df[f"{col}_found"], equal_nan=True)).sum()
        if num_mismatches != 0:
            errors.append(f'Found {num_mismatches} mismatched values while comparing column "{col}".')

    if len(errors) != 0:
        raise RuntimeError("data validation failed:\n - " + "\n - ".join(errors))

    logging.debug("%s validation was successful", path)


def extract_set_differences(ref: Set, tgt: Set, diff: Set) -> Tuple[List, List]:
    missing_vals = []
    extra_vals = []
    for v in diff:
        if v in ref:
            missing_vals.append(v)
        else:
            assert v in tgt
            extra_vals.append(v)

    return missing_vals, extra_vals


def format_error_group_mismatch(ref: Set, tgt: Set, diff: Set) -> str:
    assert len(diff) != 0
    missing_groups, extra_groups = extract_set_differences(ref, tgt, diff)

    msg = "data validation failed:"

    if len(missing_groups) != 0:
        msg += f"\n - missing group(s): {missing_groups}"

    if len(extra_groups) != 0:
        msg += f"\n - unexpected group(s): {extra_groups}"

    return msg


def format_error_dset_mismatch(ref: Set, tgt: Set, diff: Set) -> str:
    assert len(diff) != 0
    missing_groups, extra_groups = extract_set_differences(ref, tgt, diff)

    msg = "data validation failed:"

    if len(missing_groups) != 0:
        msg += f"\n - missing dataset(s): {missing_groups}"

    if len(extra_groups) != 0:
        msg += f"\n - unexpected dataset(s): {extra_groups}"

    return msg


def compare_h5_attributes(ref: h5py.Group | h5py.Dataset, tgt: h5py.Group | h5py.Dataset):
    assert ref.name == tgt.name

    logging.debug('comparing attributes for object "%s"...', ref.name)

    attr1 = dict(ref.attrs)
    attr2 = dict(tgt.attrs)
    attr1.pop("source-file", None)
    attr2.pop("source-file", None)

    diff = set(attr1) ^ set(attr2)

    if len(diff) != 0:
        missing_keys = []
        extra_keys = []
        mismatched_values = []
        for k, v in diff:
            if k in attr1:
                missing_keys.append(k)
            elif k in attr2:
                extra_keys.append(k)
            else:
                mismatched_values.append((k, attr1[k], attr2[k]))

        msg = "data validation failed:"

        if len(missing_keys) != 0:
            msg += f"\n - missing keys: {missing_keys}"

        if len(extra_keys) != 0:
            msg += f"\n - unexpected keys: {extra_keys}"

        if len(mismatched_values) != 0:
            keys = [x[0] for x in mismatched_values]
            expected_values = [x[1] for x in mismatched_values]
            found_values = [x[2] for x in mismatched_values]
            msg += f"\n - mismatched values for keys: {keys}; expected values: {expected_values}; values found: {found_values}"

        raise RuntimeError(msg)


def compare_h5_datasets(ref: h5py.Dataset, tgt: h5py.Dataset):
    assert ref.name == tgt.name
    logging.debug('comparing dataset "%s"...', ref.name)

    if ref.dtype != tgt.dtype:
        raise RuntimeError(f"dataset {ref.name} dtype is not correct: expected {ref.dtype}, found {tgt.dtype}")

    if ref.shape != tgt.shape:
        raise RuntimeError(f"dataset {ref.name} shape is not correct: expected {ref.shape}, found {tgt.shape}")

    if ref.size != tgt.size:
        raise RuntimeError(f"dataset {ref.name} size is not correct: expected {ref.size}, found {tgt.size}")

    if np.issubdtype(ref.dtype, np.inexact):
        num_mismatches = (~np.isclose(ref, tgt, equal_nan=True)).sum()
    else:
        num_mismatches = (ref[:] != tgt[:]).sum()

    if num_mismatches != 0:
        raise RuntimeError(f"dataset {ref.name} (dtype={ref.dtype}) is incorrect: found {num_mismatches} differences")

    compare_h5_attributes(ref, tgt)


def compare_h5_nodes_recursive(ref: h5py.Group, tgt: h5py.Group):
    assert ref.name == tgt.name
    logging.debug('comparing group "%s"...', ref.name)

    expected_grps = set((k for k, v in ref.items() if isinstance(v, h5py.Group)))
    found_grps = set((k for k, v in tgt.items() if isinstance(v, h5py.Group)))

    diff = expected_grps ^ found_grps
    if len(diff) != 0:
        raise RuntimeError(format_error_group_mismatch(expected_grps, found_grps, diff))

    expected_dsets = set((k for k, v in ref.items() if isinstance(v, h5py.Dataset)))
    found_dsets = set((k for k, v in tgt.items() if isinstance(v, h5py.Dataset)))

    diff = expected_dsets ^ found_dsets
    if len(diff) != 0:
        raise RuntimeError(format_error_dset_mismatch(expected_dsets, found_dsets, diff))

    compare_h5_attributes(ref, tgt)

    for dset in expected_dsets:
        compare_h5_datasets(ref[dset], tgt[dset])

    for grp in expected_grps:
        compare_h5_nodes_recursive(ref[grp], tgt[grp])


def validate_expected_profile_h5(expected: pathlib.Path, found: pathlib.Path):
    logging.debug("validating %s...", found)

    with h5py.File(expected) as ref, h5py.File(found) as tgt:
        compare_h5_attributes(ref, tgt)
        compare_h5_nodes_recursive(ref["/"], tgt["/"])


def validate_pvalue_table(expected_path: pathlib.Path, found_path: pathlib.Path):
    assert expected_path.is_file()
    try:
        if not found_path.is_file():
            raise RuntimeError(f"unable to open file {found_path}")

        expected = import_table(expected_path)
        found = import_table(found_path)

        validate_columns(expected, found, found_path)
        validate_pvalue_tables(expected, found, found_path)
    except RuntimeError as e:
        raise RuntimeError(f"failed to validate table {found_path.stem}: {e}")


def validate_nchg_cartesian_product(test_file: pathlib.Path, ref_file: pathlib.Path) -> int:
    logging.info(f"### NCHG cartesian-product: validating {test_file}...")
    if test_file.resolve() == ref_file.resolve():
        raise RuntimeError(f"test-bedpe and ref-bedpe point to the same file: {ref_file}")

    expected = import_bedpe(ref_file)
    found = import_bedpe(test_file)

    ok = True
    if (expected.columns != found.columns).any():
        logging.error("table is not correct: expected %s columns, found %s", expected.columns, found.columns)
        ok = False

    if len(expected) != len(found):
        logging.error("table is not correct: expected %d rows, found %d", len(expected), len(found))
        ok = False

    if not ok:
        return 1

    errors = []
    for col in expected.columns:
        num_mismatches = (expected[col] != found[col]).sum()
        if num_mismatches != 0:
            errors.append(f'Found {num_mismatches} mismatched values while comparing column "{col}".')

    if len(errors) != 0:
        logging.error("data validation failed:\n - " + "\n - ".join(errors))
        return 1

    logging.debug("%s validation was successful", test_file)
    return 0


def validate_nchg_compute(test_prefix: pathlib.Path, ref_prefix: pathlib.Path) -> int:
    logging.info(f"### NCHG compute: validating {test_prefix}...")
    if test_prefix == ref_prefix:
        raise RuntimeError(f"test-prefix and ref-prefix point to the same files: {ref_prefix}")

    expected_chrom_sizes = import_chrom_sizes(pathlib.Path(f"{ref_prefix}.chrom.sizes"))

    chrom_sizes_path = pathlib.Path(f"{test_prefix}.chrom.sizes")
    try:
        found_chrom_sizes = import_chrom_sizes(chrom_sizes_path)
    except Exception as e:
        logging.error(f"failed to import chrom.sizes from {chrom_sizes_path}: {e}")
        return 1

    validate_chrom_sizes(expected_chrom_sizes, found_chrom_sizes, chrom_sizes_path)

    parent_dir = ref_prefix.parent
    prefix = ref_prefix.stem

    paths = list(sorted(parent_dir.glob(f"{prefix}*.parquet")))

    if len(paths) == 0:
        raise RuntimeError(f"unable to find any tables under prefix {ref_prefix}")

    logging.debug("enumerated %d table files under prefix %s...", len(paths), ref_prefix)

    ok = True
    for path1 in paths:
        prefix = path1.as_posix().removeprefix(ref_prefix.as_posix()).removesuffix(".parquet")
        logging.debug("processing table %s...", prefix.lstrip("."))
        try:
            path2 = pathlib.Path(f"{test_prefix}{prefix}.parquet")
            validate_pvalue_table(path1, path2)
        except RuntimeError as e:
            ok = False
            logging.error(e)

        logging.debug("done processing table %s!", prefix.lstrip("."))

    return not ok


def validate_nchg_merge(test_file: pathlib.Path, ref_file: pathlib.Path) -> int:
    logging.info(f"### NCHG merge: validating {test_file}...")
    if test_file.resolve() == ref_file.resolve():
        raise RuntimeError(f"test-parquet and ref-parquet point to the same file: {ref_file}")

    try:
        validate_pvalue_table(ref_file, test_file)
        return 0
    except RuntimeError as e:
        logging.error(e)
        return 1


def validate_nchg_filter(test_file: pathlib.Path, ref_file: pathlib.Path) -> int:
    logging.info(f"### NCHG filter: validating {test_file}...")
    if test_file.resolve() == ref_file.resolve():
        raise RuntimeError(f"test-parquet and ref-parquet point to the same file: {ref_file}")

    try:
        validate_pvalue_table(ref_file, test_file)
        return 0
    except RuntimeError as e:
        logging.error(e)
        return 1


def validate_nchg_view(test_file: pathlib.Path, ref_file: pathlib.Path) -> int:
    logging.info(f"### NCHG view: validating {test_file}...")
    if test_file.resolve() == ref_file.resolve():
        raise RuntimeError(f"test-tsv and ref-tsv point to the same file: {ref_file}")

    try:
        validate_pvalue_table(ref_file, test_file)
        return 0
    except RuntimeError as e:
        logging.error(e)
        return 1


def validate_nchg_expected(test_file: pathlib.Path, ref_file: pathlib.Path) -> int:
    logging.info(f"### NCHG expected: validating {test_file}...")

    try:
        validate_expected_profile_h5(ref_file, test_file)
        return 0
    except RuntimeError as e:
        logging.error(e)
        return 1


def main() -> int:
    args = vars(make_cli().parse_args())

    setup_logger(args["verbosity"].upper())

    cmd = args["command"]
    if cmd == "cartesian-product":
        return validate_nchg_cartesian_product(args["test-bedpe"], args["ref-bedpe"])
    if cmd == "compute":
        return validate_nchg_compute(args["test-prefix"], args["ref-prefix"])

    if cmd == "merge":
        return validate_nchg_merge(args["test-parquet"], args["ref-parquet"])

    if cmd == "filter":
        return validate_nchg_filter(args["test-parquet"], args["ref-parquet"])

    if cmd == "view":
        return validate_nchg_view(args["test-tsv"], args["ref-tsv"])

    if cmd == "expected":
        return validate_nchg_expected(args["test-h5"], args["ref-h5"])

    raise NotImplementedError


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    ec = main()
    if ec == 0:
        logging.info("### status: SUCCESS!")
    else:
        logging.error("### status: FAILURE!")
    sys.exit(ec)
