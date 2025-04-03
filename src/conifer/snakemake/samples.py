"""Snakemake samples functions"""

import logging
import pathlib
import re

import pandas as pd

logger = logging.getLogger(__name__)


def samplesheet(wildcards):
    """Return samplesheet file for wildcards.sampleset"""
    resources = re.sub("data", "resources", wildcards.data)
    samplesheet_fn = pathlib.Path(
        f"{resources}/samplesets/samples-{wildcards.sampleset}.tsv"
    )
    return samplesheet_fn


def read_sampleset(wildcards, fmt="samples-{sampleset}.tsv"):
    """Read sampleset file"""
    try:
        resources = re.sub("data", "resources", wildcards.data)
    except AttributeError:
        resources = "resources"
    d = dict(wildcards)
    if "sampleset" not in d.keys():
        d["sampleset"] = "ALL"
    samplesheet_fn = pathlib.Path(f"{resources}/samplesets", fmt.format(**d))

    if not samplesheet_fn.exists():
        logger.error("No samplesheet found: %s", samplesheet_fn)
        raise FileNotFoundError
    sample_df = pd.read_csv(samplesheet_fn, sep="\t")
    return sample_df


def samplenames(wildcards):
    """Return list of sample names"""
    sample_df = read_sampleset(wildcards)
    return sample_df.SM.tolist()
