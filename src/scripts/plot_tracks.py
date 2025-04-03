#!/usr/bin/env python
"""Plot genome tracks of input BEDGRAPH files.

Plot genome tracks of input BEDGRAPH files. Input files and tracks are
specified in an input configuration file in YAML format.

The configuration file has the following format:

    reference: <path to reference FASTA index file (.fai)>
    window_size: <window size in bp>
    tracks:
      - name: <track name>
        label: <track axis label>
        file: <path to input file>
        group: <group name>
        dtype: <data type, int or float>

The tracks object is a list of dictionaries, each specifying a track
to plot. The name, label, file, and dtype keys are required. The group
key is currently unused.

"""

import argparse
import logging
import re
import sys
import textwrap
from typing import List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import yaml
from jsonschema import validate
from plotly.subplots import make_subplots

logger = logging.getLogger(__name__)

log_levels = {
    0: "WARNING",
    1: "INFO",
    2: "DEBUG",
    3: "DEBUG",
}
colors = {"0": "#1b9e77", "1": "#d95f02"}


schema = {
    "type": "object",
    "properties": {
        "reference": {
            "type": "string",
            "pattern": r"^[_/\.\-A-Za-z0-9]+.fai$",
            "description": "Path to reference FASTA index file (.fai)",
        },
        "window_size": {"type": "integer"},
        "tracks": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "label": {
                        "type": "string",
                        "description": "Track axis label",
                    },
                    "file": {"type": "string"},
                    "dtype": {"type": "string", "enum": ["int", "float"]},
                    "group": {"type": "string"},
                },
                "required": ["name", "label", "file", "dtype"],
            },
        },
    },
    "required": ["reference", "window_size", "tracks"],
    "additionalProperties": False,
}


def _sizeof_fmt(num, suffix="B"):
    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1024.0:
            return f"{num:.2f} {unit}{suffix}"
        num /= 1024.0
    return f"{num:.2f} Y{suffix}"


def _size(obj):
    return sys.getsizeof(obj)


def size(obj, humanize=True):
    """Return the size of an object"""
    if humanize:
        return _sizeof_fmt(_size(obj))
    return _size(obj)


def setup_logging(log_level) -> None:
    """Setup logging configuration"""
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def setup_parser() -> argparse.ArgumentParser:
    """Setup command line parser"""
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "config",
        metavar="CONFIG",
        help="Input configuration file in YAML format",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "output_file",
        metavar="OUTPUT_FILE",
        help="Html output file",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "--title",
        "-t",
        help="Plot title",
        type=str,
    )
    parser.add_argument(
        "--xlab",
        "-x",
        help="X-axis label",
        type=str,
    )
    parser.add_argument(
        "--ylab",
        "-y",
        help="Y-axis label",
        type=str,
    )
    # TODO: add option to specify region
    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase verbosity level. Can be supplied multiple times.",
        action="count",
        default=0,
    )

    return parser


def load_config(config_file: str) -> dict:
    """Load input configuration file"""
    logger.debug("Loading configuration file %s", config_file)
    with open(config_file) as f:  # pylint: disable=unspecified-encoding
        config = yaml.safe_load(f)

    validate(instance=config, schema=schema)

    logger.debug("Loaded configuration: %s", config)
    return config


def load_reference(reference_file: str) -> tuple:
    """Load reference FASTA index file"""
    logger.debug("Loading reference file %s", reference_file)
    reference = []
    with open(reference_file) as f:  # pylint: disable=unspecified-encoding
        for line in f:
            line = line.rstrip()
            chrom, length, offset, line_length, line_bytes = line.split("\t")
            reference.append(
                (
                    chrom,
                    {
                        "length": int(length),
                        "offset": int(offset),
                        "line_length": int(line_length),
                        "line_bytes": int(line_bytes),
                    },
                )
            )
    logger.debug("Loaded reference")
    return tuple(reference)


def get_integer_type(value) -> int:
    """Get smallest integer type that can hold value"""
    if value <= np.iinfo(np.int8).max:
        dtype = np.int8
    elif value <= np.iinfo(np.int16).max:
        dtype = np.int16
    elif value <= np.iinfo(np.int32).max:
        dtype = np.int32
    else:
        dtype = np.int64
    return dtype


def get_float_type(value) -> float:
    """Get smallest floating point type that can hold value"""
    if value <= np.finfo(np.float16).max:
        dtype = np.float16
    elif value <= np.finfo(np.float32).max:
        dtype = np.float32
    else:
        dtype = np.float64
    return dtype


def make_chromosome_levels(
    df: pd.DataFrame, unique: bool = False
) -> pd.DataFrame:
    """Make chromosome levels for coloring"""
    chrom = df.chrom
    if unique:
        levels = (chrom.cat.codes.unique() % 2).astype("str")
    else:
        levels = (chrom.cat.codes % 2).astype("str")
    return levels, chrom


def convert_to_si_suffix(number):
    """Convert a number to a string with an SI suffix."""
    suffixes = [" ", "kbp", "Mbp", "Gbp", "Tbp"]
    power = len(str(number)) // 3

    return f"{number / 1000 ** power:.1f}{suffixes[power]}"


def load_track(track: dict) -> dict:
    """Load track data from input file"""
    logger.info("Loading track %s from file %s", track["name"], track["file"])
    results = track
    names = ["chrom", "begin", "end", "value"]
    dtype = {"chrom": "category", "begin": np.int64, "end": np.int64}
    df = pd.read_csv(
        track["file"], sep="\t", header=None, names=names, dtype=dtype
    )
    vals = df["value"].describe()
    if track["dtype"] == "int":
        if vals["min"] < 0:
            raise ValueError("Negative values not allowed for integer tracks")
        if vals["max"] > 2**31 - 1:
            raise ValueError("Integer overflow")
        dtype = get_integer_type(vals["max"])
        df["value"] = df["value"].astype(dtype)
    elif track["dtype"] == "float":
        dtype = get_float_type(vals["max"])
        df["value"] = df["value"].astype(dtype)
    if "group" not in results:
        results["group"] = "default"
    results["data"] = df
    return results


def reorder_chromosomes(df: pd.DataFrame, reference: dict) -> None:
    """Reorder chromosomes in dataframe"""
    logger.info("Reordering chromosomes")
    df["chrom"] = df["chrom"].astype("category")
    desired_order = [x[0] for x in reference]
    df["chrom"] = df["chrom"].cat.reorder_categories(desired_order)


def add_track(
    fig: go.Figure,
    track: pd.DataFrame,
    row: int = 1,
) -> None:
    """Add track to figure"""
    logger.info("Adding track %s", track["name"])
    df = track["data"]
    width = df["end"] - df["begin"]
    # NB: this assumes that all windows are present in the input file
    df["genomepos"] = width.cumsum() - width.iloc[0]
    levels, chrom = make_chromosome_levels(df, unique=True)
    colormap = dict(
        zip(chrom.cat.categories.values, [colors[x] for x in levels])
    )
    line_fig = px.line(
        df,
        x="genomepos",
        y="value",
        color="chrom",
        color_discrete_map=colormap,
        labels={"genomepos": "Genome position (Gbp)", "value": track["label"]},
    )
    fig.add_trace(line_fig.data[0], row=row, col=1)


def resize_dataframe(df: pd.DataFrame) -> None:
    """Resize dataframe to reduce memory usage"""
    logger.info("Resizing data frame")
    for col in df.columns:
        if df[col].dtype == "float64":
            df[col] = df[col].astype(get_float_type(df[col].max()))
        elif df[col].dtype == "int64":
            df[col] = df[col].astype(get_integer_type(df[col].max()))


def init_track_container(reference: dict, window_size: int) -> pd.DataFrame:
    """Initialize track container"""
    logger.info("Initializing track container")
    chrom = []
    begin = []
    end = []
    for chromosome, info in reference:
        window_begin = np.arange(info["length"], step=window_size)
        begin.extend(window_begin)
        chrom.extend([chromosome] * len(window_begin))
        if len(window_begin) > 1:
            end_ar = np.concatenate(
                (window_begin[1:], np.array([info["length"]]))
            )
        elif len(window_begin) == 1:
            end_ar = np.array([info["length"]])
        else:
            end_ar = []
        if len(end_ar) != len(window_begin):
            print(len(end_ar), len(window_begin))
            raise ValueError(
                f"Length mismatch for end_ar, chromosome {chromosome}"
            )
        end.extend(end_ar)

    df = pd.DataFrame({"chrom": chrom, "begin": begin, "end": end})
    reorder_chromosomes(df, reference)

    width = df["end"] - df["begin"]
    # NB: this assumes that all windows are present in the input file
    df["genomepos"] = width.cumsum() - width.iloc[0]
    resize_dataframe(df)
    return df


def merge_tracks(
    tracks: List[pd.DataFrame], reference: dict, window_size: int
) -> pd.DataFrame:
    """Merge tracks into single dataframe"""
    logger.info("Merging tracks")
    df = init_track_container(reference, window_size)
    for track in tracks:
        logger.info("Merging track %s", track["name"])
        df = pd.merge(
            df, track["data"], how="left", on=["chrom", "begin", "end"]
        )
        df = df.rename(columns={"value": track["name"]})
    return df


def plot_tracks(
    tracks: List[pd.DataFrame],
) -> go.Figure:
    """Plot genome tracks"""
    logger.info("Plotting tracks")
    fig = make_subplots(
        rows=len(tracks), cols=1, shared_xaxes=True, vertical_spacing=0.02
    )
    for i, track in enumerate(tracks):
        logger.info("Plotting track %s", track["name"])
        add_track(fig, track, row=i + 1)
    return fig


def main(argv: List[str] = sys.argv[1:]) -> None:
    """Main entry point allowing external calls

    :param argv: command line parameter list
    :type argv: list
    """
    parser = setup_parser()
    args = parser.parse_args(argv)

    log_level = log_levels[min(args.verbose, max(log_levels.keys()))]

    setup_logging(log_level)

    logger.info("Starting plot-tracks")
    logger.debug("Arguments: %s", argv)

    config = load_config(args.config.name)

    reference = load_reference(config["reference"])

    tracks = []
    for track in config["tracks"]:
        tracks.append(load_track(track))

    data = merge_tracks(tracks, reference, config["window_size"])
    logger.info("Merged tracks into data structure of size %s", size(data))
    logger.debug("Merged tracks, head: %s", data.head())

    matrix_output = re.sub(r".html$", ".matrix.html", args.output_file.name)
    matrix_fig = px.scatter_matrix(data, dimensions=data.columns[4:])
    matrix_fig.write_html(matrix_output)

    df = (
        data.set_index(["chrom", "begin", "end", "genomepos"])
        .stack()
        .reset_index()
    )
    df.columns = ["chrom", "begin", "end", "genomepos", "plot", "value"]
    reorder_chromosomes(df, reference)
    resize_dataframe(df)
    logger.info("Stacked data frame size: %s", size(df))
    logger.debug("Stacked data summary: %s", df.info())
    logger.debug("%s", df.head())

    levels, _ = make_chromosome_levels(df, unique=True)
    colormap = dict(
        zip(df.chrom.cat.categories.values, [colors[x] for x in levels])
    )
    fig = px.scatter(
        df,
        x="genomepos",
        y="value",
        facet_row="plot",
        color="chrom",
        color_discrete_map=colormap,
    )
    fig = fig.update_yaxes(matches=None)
    args.output_file.mode = "w"
    fig.write_html(args.output_file)

    logger.info("Finished plot-tracks")


if __name__ == "__main__":
    main()
