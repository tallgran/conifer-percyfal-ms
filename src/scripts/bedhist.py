"""Summarize BEDFILE to histogram

There is a bug in d4tools in that it cannot handle disjoint BED
intervals. This script will convert a BED file to a histogram.
"""

import argparse
import logging
import sys

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
)


def main():
    """bedhist main"""
    parser = argparse.ArgumentParser(
        description="Convert BEDFILE to histogram",
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
        "--output-file",
        dest="output_file",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file name. If not provided, output is written to stdout.",
    )
    parser.add_argument(
        "--max-bins",
        dest="max_bins",
        type=int,
        default=1000,
        help="Maximum number of bins",
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    hist = np.zeros(args.max_bins + 3, dtype=int)

    while True:
        line = args.input_file.readline()
        if not line:
            break
        if line.startswith("#"):
            continue
        _, start, stop, value = line.strip().split("\t")
        start, stop = int(start), int(stop)
        try:
            hist[int(value) + 1] += stop - start
        except IndexError:
            hist[-1] += stop - start

    # Write the output data
    df = pd.DataFrame(hist)
    df.index = (
        ["<0"] + list(range(args.max_bins + 1)) + [f">{args.max_bins + 1}"]
    )

    try:
        df.to_csv(args.output_file, sep="\t", index=True, header=None)
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    main()
