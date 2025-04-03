"""Plot custom plots that don't fit in the standard plots.

This module contains custom plots that don't fit in the standard
plots. They address specific use cases (e.g., manuscript figures) and
are not intended to be used as general purpose plots.

"""

import gzip
import logging
import sys

import click
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from conifer.io.files import is_gzip

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
)


colors = {"0": "#1b9e77", "1": "#d95f02"}


def transform_genic_regions(value, size=1000):
    """Transform genic region value to integer."""
    value = value.split("_")[-1]
    fac = {"U": -1, "D": 1}[value[0]]
    i = int(value[1:])
    return fac * i * size + int(fac * size / 2)


def bootstrap(
    df,
    group_name="genic_regions",
    column="site_score",
    replicates=1000,
    confidence_interval=0.95,
):
    """Bootstrap a dataframe."""

    def _bootstrap(x):
        x = x.dropna()
        val = np.random.choice(x, size=replicates, replace=True)
        return val

    bootstrap_values = df.groupby(group_name)[column].apply(_bootstrap)
    data = []
    for group, values in bootstrap_values.items():
        ci = np.percentile(
            values,
            [
                100 * (1 - confidence_interval) / 2,
                100 * (1 + confidence_interval) / 2,
            ],
        )
        data.append([group, np.mean(values), ci[0], ci[1]])
    data = pd.DataFrame(
        data, columns=["distance", "mean", "ci_low", "ci_high"]
    )
    return data


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """Plot cli entrypoint"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__


@cli.command()
@click.argument("infile", type=click.File("rb"))
@click.option(
    "--outfile",
    "-o",
    type=click.File("wb"),
    default=sys.stdout,
    required=False,
)
@click.option("--title", type=str, default="Genic region analysis")
@click.option("--y-range", type=(float, float), default=(0.0, 0.04))
@click.option("--bootstrap-replicates", "-b", type=int, default=1000)
@click.option("--confidence-interval", "-c", type=float, default=0.95)
@click.option("--window-size", "-w", type=int, default=1000)
@click.option("--min-accessible", "-n", type=int, default=500)
def genic_region(
    infile,
    outfile,
    title,
    y_range,
    bootstrap_replicates,
    confidence_interval,
    window_size,
    min_accessible,
):
    """Custom plot to plot genic region analysis results."""
    fh = gzip.open(infile, "rt") if is_gzip(infile.name) else infile
    df = pd.read_csv(
        fh,
        dtype={
            "chrom": "category",
            "feature": "category",
            "statistic": "category",
        },
    )
    df = df[df["n_accessible"] >= min_accessible]
    df["genic_regions"] = df["region"].apply(
        transform_genic_regions, size=window_size
    )
    pi = bootstrap(
        df,
        group_name="genic_regions",
        column="site_score",
        replicates=bootstrap_replicates,
        confidence_interval=confidence_interval,
    )
    theta = bootstrap(
        df,
        group_name="genic_regions",
        column="theta",
        replicates=bootstrap_replicates,
        confidence_interval=confidence_interval,
    )
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=pi["distance"], y=pi["mean"], mode="lines", name="Mean (pi)"
        )
    )
    fig.add_trace(
        go.Scatter(
            x=pi["distance"],
            y=pi["ci_low"],
            mode="lines",
            name="CI low",
            line={"width": 0},
        )
    )
    fig.add_trace(
        go.Scatter(
            x=pi["distance"],
            y=pi["ci_high"],
            mode="lines",
            name="CI high",
            line={"width": 0},
            fill="tonexty",
            fillcolor="rgba(0,100,80,0.2)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=theta["distance"],
            y=theta["mean"],
            mode="lines",
            name="Mean (theta)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=theta["distance"],
            y=theta["ci_low"],
            mode="lines",
            name="CI low",
            line={"width": 0},
        )
    )
    fig.add_trace(
        go.Scatter(
            x=theta["distance"],
            y=theta["ci_high"],
            mode="lines",
            name="CI high",
            line={"width": 0},
            fill="tonexty",
            fillcolor="rgba(0,80,100,0.2)",
        )
    )
    fig.update_layout(
        yaxis_title="Genetic variation (%/bp)",
        xaxis_title="Distance from gene (kbp)",
        yaxis={"range": y_range},
        title=title,
    )

    outfile.mode = "w"
    fig.write_html(outfile)
