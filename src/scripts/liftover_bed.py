"""Liftover BED coordinates

There are currently two different coordinate systems in use. Due to
limitations of some variant callers to handle large chromosomes, the
chromosomes of the initial assembly was split into chunks.
This script converts between the two coordinate systems.
"""

import argparse
import gzip
import logging
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
)


def _make_coords(df):
    # NB: tacitly assumes bed file is sorted!
    start = df.iloc[0].start
    val = df.width.cumsum() + start
    x = np.append(np.array(start), val)
    return x


def _make_begin(df):
    x = _make_coords(df)
    return x[:-1]


def _make_end(df):
    x = _make_coords(df)
    return x[1:]


def _make_chromosome_index(pos: int, chunk_size: int) -> tuple:
    """Make chromosome index"""
    index = (pos // chunk_size) + 1
    new_pos = pos - (index - 1) * chunk_size
    return index, new_pos


def _process_bed6_data(data, regex, chunk_size, chunkify=False):  # pylint: disable=unused-argument
    if chunkify:
        logger.error("NOT IMPLEMENTED")
        logger.info("chunkifying data")
        return None

    logger.info("dechunkifying data")
    data["chrom"] = data.chrom_name.str.replace(regex, "", regex=True)
    data["width"] = data.stop - data.start
    data["begin"] = np.concatenate(
        data.groupby("chrom").apply(_make_begin).values
    )
    data["end"] = np.concatenate(data.groupby("chrom").apply(_make_end).values)
    return data[["chrom", "begin", "end", "feature", "strand", "value"]]


def _process_bed_data(data, regex, chunk_size, chunkify=False):  # pylint: disable=unused-argument
    if chunkify:
        logger.error("NOT IMPLEMENTED")
        logger.info("chunkifying data")
        return None
    logger.info("dechunkifying data")
    data["chrom"] = data.chrom_name.str.replace(regex, "", regex=True)
    data["width"] = data.stop - data.start
    data["begin"] = np.concatenate(
        data.groupby("chrom").apply(_make_begin).values
    )
    data["end"] = np.concatenate(data.groupby("chrom").apply(_make_end).values)
    return data[["chrom", "begin", "end", "value"]]


def _process_gff_data(data, regex, chunk_size, chunkify=False):
    if chunkify:
        logger.info("chunkifying data")
        x = data.start.apply(_make_chromosome_index, chunk_size=chunk_size)
        data["index"] = [f"_{i[0]}" for i in x]
        new_pos = [i[1] for i in x]
        index_count = defaultdict(list)
        for key, value in zip(data["chrom_name"], data["index"]):
            if value not in index_count[key]:
                index_count[key].append(value)
        zero_index = [
            key for key, value in index_count.items() if len(value) == 1
        ]
        subset = data.loc[data["chrom_name"].isin(zero_index)].copy()
        subset["index"] = ""
        data.update(subset)
        data["chrom"] = data["chrom_name"].str.cat(data["index"].astype(str))
        data["start"] = new_pos
        # End may be in a different chunk
        x = data.stop.apply(_make_chromosome_index, chunk_size=chunk_size)
        new_pos = [i[1] for i in x]
        data["stop"] = np.array(new_pos).astype(int)
        subset = data.iloc[np.where(data["stop"] < data["start"])].copy()
        subset["stop"] = chunk_size
        data.update(subset)
        data["start"] = data["start"].astype(int)
        data["stop"] = data["stop"].astype(int)
    else:
        logger.info("dechunkifying data")
        data["chrom"] = data.chrom_name.str.replace(regex, "", regex=True)
        data["width"] = data.stop - data.start
        data["start"] = np.concatenate(
            data.groupby("chrom").apply(_make_begin).values
        )
        data["stop"] = np.concatenate(
            data.groupby("chrom").apply(_make_end).values
        )
    return data[
        [
            "chrom",
            "source",
            "feature",
            "start",
            "stop",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
    ]


def process_data(
    data: pd.DataFrame, regex: str, chunk_size: int, chunkify: bool = False
) -> pd.DataFrame:
    """Process BED/GFF data"""
    if data.shape[1] == 9:
        logger.info("Processing GFF data")
        return _process_gff_data(data, regex, chunk_size, chunkify)
    if data.shape[1] == 4:
        logger.info("Processing BED data")
        return _process_bed_data(data, regex, chunk_size, chunkify)
    if data.shape[1] == 6:
        logger.info("Processing BED6 data")
        return _process_bed6_data(data, regex, chunk_size, chunkify)

    logger.error("Data format not recognized")
    return None


def main():
    """liftoverBED main function"""
    parser = argparse.ArgumentParser(
        description="Chunkify/dechunkify BED format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add an argument for file input
    parser.add_argument(
        "input_file",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input file name. If not provided, input is read from stdin.",
    )

    # Add an argument for file output
    parser.add_argument(
        "-o",
        "--output_file",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help=(
            "Output file name. If "
            "not provided, output is written to stdout."
        ),
    )
    parser.add_argument(
        "--chunk-size", type=int, default=100000000, help="Chunk size"
    )
    parser.add_argument(
        "--chunkify",
        action="store_true",
        help=("Chunkify input. Will add UNDERSCORE_DIGIT to chromosome names"),
    )
    parser.add_argument(
        "--regex",
        type=str,
        default=r"_\d+$",
        help=(
            "regex that identifies chunks. This will "
            "be stripped from chromosome names"
        ),
    )
    parser.add_argument(
        "--type",
        type=str,
        default="bed",
        choices=["bed", "gff"],
        help="Input file type. Options: bed, gff",
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Read the input data
    if args.type == "bed":
        colnames = ["chrom_name", "start", "stop", "value"]
    else:
        colnames = [
            "chrom_name",
            "source",
            "feature",
            "start",
            "stop",
            "score",
            "strand",
            "frame",
            "attribute",
        ]

    if args.input_file.name.endswith(".gz"):
        df = pd.read_table(
            gzip.open(args.input_file.name, "rt"),
            header=None,
            comment="#",
            names=colnames,
        )
    else:
        df = pd.read_table(
            args.input_file,
            header=None,
            comment="#",
            names=colnames,
        )

    # Process input
    df = process_data(
        df, rf"{args.regex}", args.chunk_size, chunkify=args.chunkify
    )

    # Write the output data
    try:
        df.to_csv(args.output_file, sep="\t", index=False, header=None)
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    main()
