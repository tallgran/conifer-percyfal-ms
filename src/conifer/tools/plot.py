"""Plotting utilities"""

import gzip
import logging
import os
import pathlib
import re
import sys

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import seaborn as sns
from plotly.subplots import make_subplots

from conifer.io.files import is_gzip, sniff_infile

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)

log_levels = {
    0: "WARNING",
    1: "INFO",
    2: "DEBUG",
    3: "DEBUG",
}
colors = {"0": "#1b9e77", "1": "#d95f02"}


def setup_logging(log_level=1) -> None:
    """Setup logging configuration"""
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_bedgraph(file):
    """Read bedgraph into dataframe"""
    delimiter, _ = sniff_infile(file.name)
    data = pd.read_csv(
        file.name,
        header=None,
        names=["chrom_name", "start", "stop", "value"],
        sep=delimiter,
    )
    data["chrom"] = data["chrom_name"].apply(lambda x: re.sub(r"_\d+$", "", x))
    data["width"] = data.stop - data.start

    def _make_pos(df):
        wzero = df.iloc[0].width
        start = df.iloc[0].start - wzero
        val = df.width.cumsum() + start
        return val

    data["pos"] = data.groupby("chrom").apply(_make_pos).values
    data["genomepos"] = data.pos.cumsum()
    return data


def read_hist_data(fn, label, show_feature_size=True, max_bin=None):
    """Read histogram data"""
    logger.info("Reading %s", fn.name)
    data = pd.read_table(fn, header=None, names=["coverage", "count"])
    data["fn"] = [fn.name] * data.shape[0]
    data["x"] = data["coverage"].str.replace(">", "")
    data.at[0, "x"] = "-1"
    data["x"] = data["x"].astype(int)
    data["count"] = data["count"].astype(int)
    if max_bin is not None:
        if max_bin > 1:
            max_bin = int(max_bin)
        else:
            pct = np.cumsum(data["count"]) / max(np.cumsum(data["count"]))
            max_bin = np.min(np.where(pct > max_bin))
        data.at[max_bin, "count"] = np.sum(data["count"].values[max_bin:])
        data = data.drop(range((max_bin + 1), data.shape[0]))
        cov = data["coverage"].values[-1]
        data.at[data.shape[0] - 1, "coverage"] = f">{cov}"
    if show_feature_size:
        label = f"{label} ({convert_to_si_suffix(np.sum(data['count']))})"
    data["label"] = [label] * data.shape[0]
    return data


def make_chromosome_levels(data, unique=False):
    """Make chromosome levels for coloring"""
    chrom = data.chrom.astype("category")
    if unique:
        levels = (chrom.cat.codes.unique() % 2).astype("str")
    else:
        levels = (chrom.cat.codes % 2).astype("str")
    return levels, chrom


def convert_to_si_suffix(number):
    """Convert a number to a string with an SI suffix."""
    suffixes = [" ", "kbp", "Mbp", "Gbp", "Tbp"]
    power = len(str(int(number))) // 3
    return f"{number / 1000 ** power:.1f}{suffixes[power]}"


def check_outfile(outfile):
    """Check that outfile does not exist"""
    if pathlib.Path(outfile).exists():
        logger.error("%s exists!", outfile)
        logger.error("Make sure to provide a non-existing output file name")
        sys.exit()


def make_vector(df, sample_size=10000):
    """Make subsampled vector from dataframe"""
    n = np.sum(df["count"])
    return np.random.choice(
        df["x"], size=min(sample_size, n), p=df["count"] / n
    )


def mpl_plot(data, outfile, *, title, width, height, **kwargs):  # pylint: disable=too-many-locals
    """Plot with matplotlib"""
    fig, ax = plt.subplots(figsize=(width, height))
    fig.set_size_inches(width, height)
    kind = kwargs.pop("kind", "line")
    x = kwargs.pop("x", "x")
    y = kwargs.pop("y", "y")
    if kind in ["boxplot", "violin"]:
        by = kwargs.get("by", "feature")
        column = kwargs.get("column", "count")
        features = data[by].unique()
        data = [
            data[data[by] == feature][column].values for feature in features
        ]
        if kind == "boxplot":
            plt.boxplot(data)
        elif kind == "violin":
            plt.violinplot(data)
        plt.xticks(range(1, len(features) + 1), features)
        plt.xticks(rotation=45, ha="right")
    else:
        if isinstance(data, pd.DataFrame):
            data.plot(ax=ax, **kwargs)
        else:
            num_bins = kwargs.pop("num_bins", None)
            for group, group_data in data:
                logger.info("Plotting %s", group)
                kw = {"label": group}
                if "level" in group_data.columns:
                    kw["color"] = colors.get(
                        group_data.level.unique()[0], "black"
                    )
                # Plot the data points with alternating colors
                ax.plot(group_data[x], group_data[y], **kw)
                ax.set_xlabel(kwargs.get("xlab"))
                ax.set_ylabel(kwargs.get("ylab"))
                if kwargs.get("legend", True):
                    plt.legend()
                if num_bins is not None:
                    ax2 = ax.twinx()
                    ax2.hist(
                        group_data[x],
                        weights=group_data[y],
                        bins=num_bins,
                        label=group,
                        alpha=0.2,
                        color="red",
                    )
                    ax2.set_ylabel(kwargs.get("ylab"))
    plt.title(title)

    with outfile.open() as fh:
        fig.savefig(fh, bbox_inches="tight", format="png")


def seaborn_plot(data, outfile, *, title, width, height, **kwargs):
    """Make seaborn plot"""
    fig, _ = plt.subplots(figsize=(width, height))
    fig.set_size_inches(width, height)
    sns.set_style("whitegrid")
    if "colors" in data.columns:
        sns.set_palette(data.colors.unique())
    kind = kwargs.pop("kind", "line")
    x = kwargs.pop("x", "x")
    y = kwargs.pop("y", "y")
    if kind in ["boxplot", "violin"]:
        xlab = kwargs.pop("xlab", "Feature")
        ylab = kwargs.pop("ylab", "Count")
        x = kwargs.get("by", "feature")
        y = kwargs.get("column", "count")
        if kind == "violin":
            sns.violinplot(data=data, x=x, y=y, **kwargs)
        elif kind == "boxplot":
            sns.boxplot(data=data, x=x, y=y, **kwargs)
        plt.xticks(rotation=45, ha="right")
    elif kind == "line":
        num_bins = kwargs.pop("num_bins", None)
        xlab = kwargs.pop("xlab", "x")
        ylab = kwargs.pop("ylab", "y")
        sns.relplot(
            data=data,
            x=x,
            y=y,
            height=height,
            aspect=width / height,
            kind=kind,
            **kwargs,
        )
        if num_bins is not None:
            ax2 = plt.twinx()
            sns.histplot(
                data=data,
                x=x,
                y=y,
                ax=ax2,
                bins=num_bins,
                color="red",
                alpha=0.2,
            )
            ax2.set_ylabel(ylab)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    with outfile.open() as fh:
        plt.savefig(fh, bbox_inches="tight", format="png")


def plotly_plot(data, outfile, *, title, xlab, ylab, **kwargs):
    """Make plotly plot"""
    kind = kwargs.pop("kind", "line")
    x = kwargs.pop("x", "x")
    y = kwargs.pop("y", "y")
    if kind == "line":
        num_bins = kwargs.pop("num_bins", None)
        if isinstance(data, pd.DataFrame):
            fig = px.line(
                data,
                x=x,
                y=y,
                title=title,
                **kwargs,
            )
        else:
            if num_bins is not None:
                fig = make_subplots(specs=[[{"secondary_y": True}]])
            else:
                fig = go.Figure()
            for group, group_data in data:
                logger.info("Plotting %s", group)
                fig.add_trace(
                    go.Scatter(x=group_data[x], y=group_data[y], name=group)
                )
                if num_bins is not None:
                    fig.add_trace(
                        go.Histogram(
                            x=group_data[x],
                            y=group_data[y],
                            histfunc="sum",
                            nbinsx=num_bins,
                            opacity=0.5,
                            name=group,
                        ),
                        secondary_y=True,
                    )
                    fig.update_layout(
                        yaxis2={"tickformat": ".2s"},
                    )
    elif kind in ["boxplot", "violin"]:
        fig = go.Figure()
        for group, group_data in data:
            v = make_vector(group_data)
            if kind == "boxplot":
                fig.add_trace(go.Box(y=v, name=group))
            else:
                fig.add_trace(go.Violin(y=v, name=group))
    fig.update_layout(
        title=title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        xaxis={"tickformat": ".2s"},
        yaxis={"tickformat": ".2s"},
    )
    outfile.mode = "w"
    fig.write_html(outfile)


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """Plot cli entrypoint"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__


@cli.command()
@click.argument("infile", type=click.File("r"))
@click.argument(
    "outfile", type=click.File("wb"), default=sys.stdout, required=False
)
@click.option(
    "--backend", type=click.Choice(["plotly", "mpl", "seaborn"]), default="mpl"
)
@click.option("--title", type=str)
@click.option("--width", type=float, default=9.0)
@click.option("--height", type=float, default=4.5)
@click.pass_context
def coverage(  # noqa: A001, pylint: disable=too-many-arguments,too-many-locals
    ctx,
    infile,
    outfile,
    backend,
    title,
    width,
    height,
):
    """Plot coverage over contigs"""
    setup_logging()
    logger.info(
        "Running %s version %s", ctx.find_root().info_name, ctx.obj["VERSION"]
    )
    data = read_bedgraph(infile)
    xlab = "Genome position (bp)"
    ylab = "Coverage (X)"
    if backend == "plotly":
        levels, chrom = make_chromosome_levels(data, unique=True)
        colormap = dict(
            zip(chrom.cat.categories.values, [colors[x] for x in levels])
        )
        plotly_plot(
            data,
            outfile,
            title=title,
            x="genomepos",
            xlab=xlab,
            ylab=ylab,
            y="value",
            color="chrom",
            labels={
                "genomepos": "Genome position (Gbp)",
                "value": ylab,
            },
            color_discrete_map=colormap,
        )
    elif backend == "mpl":
        levels, _ = make_chromosome_levels(data)
        data["level"] = levels
        data["colors"] = levels.apply(lambda x: colors[x])
        mpl_plot(
            data.groupby("chrom"),
            outfile,
            title=title,
            xlab=xlab,
            ylab=ylab,
            width=width,
            height=height,
            x="genomepos",
            y="value",
            color="chrom",
        )
    elif backend == "seaborn":
        levels, _ = make_chromosome_levels(data)
        data["level"] = levels
        data["colors"] = levels.apply(lambda x: colors[x])
        seaborn_plot(
            data,
            outfile,
            title=title,
            width=width,
            height=height,
            xlab="Genome position (bp)",
            ylab="Coverage (X)",
            x="genomepos",
            y="value",
            hue="chrom",
        )


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
    "--backend", type=click.Choice(["plotly", "mpl", "seaborn"]), default="mpl"
)
@click.option("--width", type=float, default=9.0, help="plot width (inches)")
@click.option("--height", type=float, default=4.5, help="plot height (inches)")
@click.option("--title", type=str, help="plot title")
@click.option("--xlab", type=str, help="x label")
@click.option("--ylab", type=str, help="y label")
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
    "--labels", type=str, default=None, help="Comma-separated list of labels"
)
@click.option(
    "--show-feature-size",
    is_flag=True,
    default=False,
    help="Show feature size",
)
@click.pass_context
def hist(  # noqa: A001, pylint: disable=too-many-arguments,too-many-locals
    ctx,
    infile,
    outfile,
    backend,
    title,
    width,
    height,
    xlab,
    ylab,
    max_bin,
    num_bins,
    labels,
    show_feature_size,
):
    """Plot histogram of coverage.

    Read in a mosdepth hist file with two columns: coverage and count
    and plot histogram as a line plot or binned bar plot.
    """
    setup_logging()
    logger.info(
        "Running %s version %s", ctx.find_root().info_name, ctx.obj["VERSION"]
    )

    if labels is not None:
        labels = labels.split(",")
    else:
        labels = [os.path.basename(fn.name) for fn in infile]
    title = title if title else "Coverage histogram"
    xlab = xlab if xlab else "Coverage"
    ylab = ylab if ylab else "Count"

    ext = os.path.splitext(outfile.name)[1].lstrip(".")
    if ext == "html":
        backend = "plotly"

    dflist = []
    for i, fn in enumerate(infile):
        dflist.append(
            read_hist_data(fn, labels[i], show_feature_size, max_bin)
        )

    data = pd.concat(dflist)
    groups = data.groupby("label")

    if backend == "plotly":
        plotly_plot(
            groups,
            outfile,
            title=title,
            x="x",
            y="count",
            xlab=xlab,
            ylab=ylab,
            num_bins=num_bins,
        )
    elif backend == "mpl":
        mpl_plot(
            groups,
            outfile,
            title=title,
            xlab=xlab,
            ylab=ylab,
            width=width,
            height=height,
            num_bins=num_bins,
            y="count",
        )
    elif backend == "seaborn":
        seaborn_plot(
            data,
            outfile,
            title=title,
            xlab=xlab,
            ylab=ylab,
            width=width,
            height=height,
            num_bins=num_bins,
            y="count",
        )


@cli.command()
@click.argument("infile", type=click.File("r"), nargs=-1)
@click.option(
    "--outfile",
    "-o",
    type=click.File("wb"),
    default=sys.stdout,
    required=False,
)
@click.option(
    "--backend", type=click.Choice(["plotly", "mpl", "seaborn"]), default="mpl"
)
@click.option("--title", type=str)
@click.option("--width", type=float, default=9.0)
@click.option("--height", type=float, default=4.5)
@click.option("--xlab", type=str, default="Feature")
@click.option("--ylab", type=str, default="Count")
@click.option("--sample-size", type=int, default=10000)
@click.option("--labels", type=str, default=None)
@click.option(
    "--show-feature-size",
    is_flag=True,
    default=False,
    help="Show feature size",
)
@click.option(
    "--plot-type", type=click.Choice(["boxplot", "violin"]), default="boxplot"
)
@click.pass_context
def boxplot_hist(  # noqa: A001, pylint: disable=too-many-arguments,too-many-locals
    ctx,
    infile,
    outfile,
    backend,
    title,
    width,
    height,
    xlab,
    ylab,
    sample_size,
    labels,
    show_feature_size,
    plot_type,
):
    """Plot boxplot of histogram files"""
    setup_logging()
    logger.info(
        "Running %s version %s", ctx.find_root().info_name, ctx.obj["VERSION"]
    )

    title = (
        f"{title} Sample size={sample_size}"
        if title
        else f"Sample size={sample_size}"
    )
    xlab = xlab if xlab else "Coverage"
    ylab = ylab if ylab else "Count"
    if labels is not None:
        labels = labels.split(",")
    else:
        labels = [os.path.basename(fn.name) for fn in infile]

    dflist = []
    for i, fn in enumerate(infile):
        dflist.append(read_hist_data(fn, labels[i], show_feature_size))

    df = pd.concat(dflist)

    groups = df.groupby("label")

    if backend == "plotly":
        plotly_plot(
            groups,
            outfile,
            title=title,
            x="x",
            y="count",
            xlab=xlab,
            ylab=ylab,
            kind=plot_type,
        )
    else:
        dflist = []
        for group, group_data in groups:
            v = make_vector(group_data)
            dflist.append(pd.DataFrame({"feature": group, "count": v}))
        df = pd.concat(dflist)
        if backend == "mpl":
            mpl_plot(
                df,
                outfile,
                title=title,
                xlab=xlab,
                ylab=ylab,
                width=width,
                height=height,
                kind=plot_type,
            )
        elif backend == "seaborn":
            seaborn_plot(
                df,
                outfile,
                title=title,
                width=width,
                height=height,
                xlab=xlab,
                ylab=ylab,
                kind=plot_type,
                hue="feature",
            )


@cli.command()
@click.argument("infile", type=click.File("rb"))
@click.option("--outfile", "-o", type=click.File("wb"))
@click.option("--title", type=str)
@click.option(
    "--ordering",
    type=str,
    default=None,
    help="Subset and ordering of features",
)
@click.option(
    "--feature-size-column",
    type=str,
    default=None,
    help="Show feature size using feature size column to compute size",
)
@click.option("--width", type=float, default=900)
@click.option("--height", type=float, default=450)
@click.option("-x", type=str, default="x", help="x column")
@click.option("-y", type=str, default="y", help="y column (features)")
@click.option(
    "--verbose", "-v", count=True, default=0, help="Increase verbosity level"
)
@click.pass_context
def boxplot_long(  # pylint: disable=too-many-arguments,too-many-locals
    ctx,
    infile,
    outfile,
    title,
    ordering,
    feature_size_column,
    x,
    y,
    verbose,
    width,  # pylint: disable=unused-argument
    height,  # pylint: disable=unused-argument
):
    """Make boxplot of long form data"""
    log_level = log_levels[min(verbose, max(log_levels.keys()))]
    setup_logging(log_level)
    logger.info(
        "Running %s version %s", ctx.find_root().info_name, ctx.obj["VERSION"]
    )
    title = title if title else "Boxplot"

    fh = gzip.open(infile, "rt") if is_gzip(infile.name) else infile
    df = pd.read_csv(fh)
    df[x] = df[x].astype("category")
    if ordering is not None:
        ordering = ordering.split(",")
        df = df[df[x].isin(ordering)]
        df[x] = df[x].cat.remove_unused_categories()
        df[x] = df[x].cat.reorder_categories(ordering, ordered=True)
    else:
        ordering = df[x].unique()

    fig = px.box(
        df,
        x=x,
        y=y,
        title=title,
        template="plotly_white",
        category_orders={x: ordering},
    )
    if feature_size_column is not None:
        feature_size = df.groupby(x)[feature_size_column].sum()
        labels = [
            f"{x} ({convert_to_si_suffix(y)})"
            for x, y in zip(df[x].cat.categories, feature_size)
        ]
        fig.update_xaxes(
            tickangle=45,
            tickvals=df[x].cat.categories,
            ticktext=labels,
            tickmode="array",
        )
    if outfile is None:
        sys.stdout.buffer.write(fig.to_image(format="svg"))
    else:
        if outfile.name.endswith(".html"):
            outfile.mode = "w"
            fig.write_html(outfile)
        elif outfile.name.endswith(".png"):
            pio.write_image(fig, outfile.name)
        else:
            logger.error("Unknown output file type %s", outfile.name)
