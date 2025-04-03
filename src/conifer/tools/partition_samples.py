r"""Partition SAMPLEFILE into smaller sample sets.

Partition SAMPLEFILE into subsets, possibly based on columns in sample
sheet, or select samples from another input file.

"""

import logging
import os
import pathlib
import re

import click
import numpy as np
import pandas as pd
from tqdm import tqdm

from .. import __version__

logger = logging.getLogger(__name__)


def partition_by_column(
    samples,
    columns,
    single,
    output_directory,
    outfile,
    partition_size,
    which,
    label,
):
    for g, x in samples.groupby(list(columns)):
        if isinstance(g, str):
            g = [g]
        elif isinstance(g, tuple):
            g = list(g)
        gg = "-".join(g).replace(" ", "-")
        n = x.shape[0]
        if single:
            if outfile is None:
                outfile = pathlib.Path(output_directory) / f"{label}-{gg}.tsv"
            x.to_csv(outfile, sep="\t", index=False)
        else:
            indices = np.arange(0, n, partition_size)
            indices = np.append(indices, n)
            for i, (j, k) in enumerate(zip(indices[0:-1], indices[1:])):
                outfile = (
                    pathlib.Path(output_directory)
                    / f"{label}-{gg}-{i+1}-of-{len(indices[0:-1])}.tsv"
                )
                if which is not None:
                    if which == i + 1:
                        x.iloc[j:k].to_csv(outfile, sep="\t", index=False)
                else:
                    x.iloc[j:k].to_csv(outfile, sep="\t", index=False)


def partition_by_coverage(
    samples, coverage_files, output_directory, outfile, covmin, covmax, label
):
    pdlist = []
    if len(coverage_files) > 0:
        for cf in tqdm(coverage_files):
            x = pd.read_table(cf)
            x["SM"] = re.sub(
                ".mosdepth.summary.txt$", "", os.path.basename(cf)
            )
            pdlist.append(x.loc[x.chrom == "total"])
        df = pd.concat(pdlist)
    x = samples.merge(df, left_on="SM", right_on="SM")

    if outfile is None:
        outfile = (
            pathlib.Path(output_directory)
            / f"{label}-coverage-min{covmin}-max{covmax}.tsv"
        )
    x.loc[(x["mean"] >= covmin) & (x["mean"] <= covmax)].to_csv(
        outfile, sep="\t", index=False
    )


@click.command(help=__doc__)
@click.version_option(version=__version__)
@click.argument("samplesheet", type=click.Path(exists=True))
@click.option(
    "--partition-size",
    "-p",
    help="requested maximum partition sizes",
    default=100,
    type=int,
    show_default=True,
)
@click.option(
    "--columns",
    "-c",
    help="columns to partition on",
    type=str,
    multiple=True,
)
@click.option(
    "--output-directory",
    "-o",
    help="output directory",
    default=pathlib.Path(os.curdir),
    type=click.Path(exists=True),
)
@click.option(
    "outfile",
    "--output-file-name",
    "-O",
    help="output file name",
    type=click.Path(exists=False),
)
@click.option(
    "--exclude-expr",
    "-e",
    help=(
        "pandas expression to exclude samples from groups "
        "expressed as column:value"
    ),
    default=None,
    type=str,
    multiple=True,
)
@click.option(
    "--include-expr",
    "-i",
    help=(
        "pandas expression to include samples from groups "
        "expressed as column:value"
    ),
    default=None,
    type=str,
    multiple=True,
)
@click.option(
    "--which",
    "-w",
    help="output just one subset; required to avoid snakemake checkpoints",
    default=None,
    type=int,
)
@click.option(
    "--single",
    "-s",
    is_flag=True,
    help=(
        "output all samples that pass criteria to one file, regardless of "
        "partition size. Output file names will not contain i-of-n string."
    ),
)
@click.option(
    "--covmin",
    help="minimum coverage",
    type=float,
)
@click.option(
    "--covmax",
    help="maximum coverage",
    type=float,
)
@click.argument(
    "coverage-files",
    type=click.Path(exists=True),
    nargs=-1,
)
@click.option(
    "--omit-samples",
    type=click.Path(exists=True),
    help="list of samples to omit from analysis",
)
def cli(
    samplesheet,
    partition_size,
    columns,
    output_directory,
    outfile,
    exclude_expr,
    include_expr,
    which,
    single,
    covmin,
    covmax,
    coverage_files,
    omit_samples,
):
    samples = pd.read_table(samplesheet)
    assert set(columns).issubset(
        samples.columns.values
    ), f"some columns not found in sample sheet: {columns}"
    label = "samples"
    if exclude_expr is not None:
        for e in exclude_expr:
            col, val = e.split(":")
            samples = samples[samples[col] != val]
    if include_expr is not None:
        for i in include_expr:
            col, val = i.split(":")
            samples = samples[samples[col] == val]
    if omit_samples is not None:
        logger.info(f"omitting samples in {omit_samples}")
        omit_samples_list = pd.read_table(
            omit_samples, comment="#", header=None
        )
        samples = samples[~samples.SM.isin(omit_samples_list[0])]
        label = (
            f"{label}-{os.path.splitext(os.path.basename(omit_samples))[0]}"
        )
    if len(columns) > 0:
        logger.info("partitioning on column(s) {columns}")
        partition_by_column(
            samples,
            columns,
            single,
            output_directory,
            outfile,
            partition_size,
            which,
            label,
        )
    if covmin is not None and covmax is not None:
        logger.info(
            f"subsetting on coverage range (inclusive): {covmin} - {covmax}"
        )
        partition_by_coverage(
            samples,
            coverage_files,
            output_directory,
            outfile,
            covmin,
            covmax,
            label,
        )
