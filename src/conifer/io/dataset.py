"""Sgkit genome dataset functions"""

import gzip
import logging
import re

import numpy as np
import pandas as pd
import xarray as xr

logger = logging.getLogger(__name__)


def create_dataset_from_fai(fai: str) -> xr.Dataset:
    """Create a genome dataset from a fasta index

    :param str fai: fasta index
    :returns: Genome dataset
    :rtype: xr.Dataset
    """
    data = pd.read_table(
        fai, names=["chrom", "length", "offset", "linebases", "linewidth"]
    )
    ds = xr.Dataset(data_vars=dict.fromkeys(data.chrom, []))
    for chrom, length in zip(data.chrom, data.length):
        ds[chrom] = np.zeros(length, dtype="<u4")
    return ds


def pivot_coverage(data):
    """Pivot coverage data"""
    x = np.zeros((max(data.chromEnd) - min(data.chromStart) + 1,), dtype="<i4")
    for start, end, cov in zip(data.chromStart, data.chromEnd, data.coverage):
        s = start + 1
        x[s:end] = cov
    return x


def get_filetype(filename):
    """Get filetype from filename"""
    if re.search(r".d4$", filename):
        return "d4"
    if re.search(r".bed(.gz)?$", filename):
        return "bed"
    logger.error("Unsupported filetype: {filename}")
    return None


def txt(filename):
    """OBSOLETE: Read coverage file"""
    cov = pd.read_table(
        filename, names=["chrom", "chromStart", "chromEnd", "coverage"]
    )
    with gzip.open("coverage.txt.gz", "wb") as fh:
        for chrom, data in cov.groupby("chrom"):
            x = pivot_coverage(data)
            fh.write(f"{chrom}\t".encode())
            np.savetxt(fh, x, newline="\t", fmt="%i")
            fh.write(b"\n")
