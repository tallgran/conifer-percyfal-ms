"""Summarize stat.

Summarize statistics from a file.
"""

import logging
import os
import sys

import click
import numpy as np
import pandas as pd
import scipy

from .. import __version__
from .plot import read_hist_data

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """summarize docstring"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.info("Running %s version %s", __name__, __version__)


@cli.command()
@click.argument("infile", type=click.File("r"), nargs=-1)
@click.option(
    "--outfile",
    "-o",
    type=click.File("wb"),
    default=sys.stdout,
    required=False,
    help="Output file name",
)
@click.option(
    "--max-bin",
    type=float,
    help=(
        "Set the upper bound of max bin. All bins above max bin "
        "will be summed. "
        "A float in the range 0-1 will be treated as a percentile."
    ),
)
@click.option(
    "--num-bins",
    type=int,
    help="Overlay a binned histogram with num_bins bins.",
)
@click.option(
    "--show-feature-size",
    is_flag=True,
    default=False,
    help="Show feature size",
)
@click.pass_context
def hist(ctx, infile, outfile, max_bin, num_bins, show_feature_size):
    """Calculate summary statistics from d4 hist file"""
    logger.info("Running hist")
    logger.info("infile: %s", infile)
    logger.info("outfile: %s", outfile)

    labels = [os.path.basename(fn.name) for fn in infile]
    dflist = []
    for i, fn in enumerate(infile):
        dflist.append(
            read_hist_data(fn, labels[i], show_feature_size, max_bin)
        )

    data = pd.concat(dflist)
    # groups = data.groupby("label")
    print(data.tail())

    p = data["count"] / np.sum(data["count"])
    x = np.random.choice(data["x"], p=p, size=10000, replace=True).astype(
        np.int64
    )
    a = np.arange(0, 1.1, 0.1)
    x = data["x"]
    print(x)
    peaks = scipy.signal.find_peaks(x, height=0)[0]
    print(np.mean(x), np.median(x), np.std(x), np.quantile(x, a))
    print(peaks)
