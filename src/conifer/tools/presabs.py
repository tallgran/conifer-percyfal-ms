"""presabs

Presence-absence analysis of coverage data.
"""

import logging

import click
import numpy as np
import pandas as pd

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """Compute presence/absence of coverage data."""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.debug("Running %s version %s", __shortname__, __version__)


@cli.command()
@click.argument("bed", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option("--threshold", "-t", help="coverage threshold", default=3)
def summarize(bed, threshold, output):
    """Summarize coverage data.

    The input BED file should consist of five columns (strand is
    excluded).
    """
    logger.info("Summarizing coverage data")
    df = pd.read_table(bed, header=None)
    if df.shape[1] != 5:
        raise ValueError("BED file must have five columns")
    df.columns = ["chrom", "start", "end", "name", "score"]
    grouped = df.groupby("name")
    with open(output, "w") as outfile:
        for name, group in grouped:
            presabs = (group["score"] > threshold).astype(np.int8)
            widths = group["end"] - group["start"]
            cov = presabs * widths
            outfile.write(f"{name}\t{cov.sum()}\t{widths.sum()}\n")
