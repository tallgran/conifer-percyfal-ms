"""Snakemake parameter functions"""

import math
import re
from pathlib import Path

from conifer.snakemake import samples


def d4_hist_maxbin(wildcards):
    """Return maxbin parameter for rule d4_hist"""
    sample_df = samples.read_sampleset(wildcards)
    if wildcards.prefix.startswith("sum"):
        maxbin = 50 * len(sample_df.SM)
        n = 10 ** (len(str(maxbin)) - 1)
        maxbin = math.ceil(maxbin / n) * n
    else:
        maxbin = len(sample_df.SM) + 1
    return maxbin


def conifer_plot_hist_num_bins(wildcards):
    """Return num_bins for rule conifer_plot_hist"""
    sample_df = samples.read_sampleset(wildcards)
    if wildcards.prefix.startswith("sum"):
        num_bins = min(20, len(sample_df.SM))
    else:
        num_bins = 50
    return num_bins


def aggregate_coverage_plot_title(wildcards, output):
    """Return plot title for aggregate coverage plots"""
    sample_df = samples.read_sampleset(wildcards)
    n = sample_df.shape[0]
    title = f"{wildcards.sampleset} (n={n}). "
    if "feature" in wildcards.keys():
        title += f"Feature: {wildcards.feature}. "
    if "prefix" in wildcards.keys():
        if wildcards.prefix.startswith("sum"):
            title += "Sum coverage."
        else:
            coverage = re.sub(
                ".p[0-9]+", "", re.sub("count_ge", "", wildcards.prefix)
            )
            title += f"Count coverage >= {coverage}X."
        return title
    outfile = Path(output.pop())
    if outfile.name.startswith("sum"):
        title += "Sum coverage"
    else:
        title += f"Count coverage >= {wildcards.coverage}X"
        try:
            title += f", window size {wildcards.window_size}."
        except AttributeError:
            title += "."
    return title
