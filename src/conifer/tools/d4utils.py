"""d4utils - utilities to work with d4 file format."""

import logging
import pathlib
import re
import sys

import click
import numpy as np
import pandas as pd
import pyd4
from numpy.lib.stride_tricks import sliding_window_view
from tqdm import tqdm

from .. import __version__

__shortname__ = __name__.split(".")[-1]

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
)

pat = re.compile(r"[ \t]+")


def make_chunks(begin, end, size):
    pos = np.arange(begin, end, size)
    begin_list = pos
    end_list = pos[1 : len(pos)]  # noqa: E203
    end_list = np.append(end_list, end)
    for begin, end in zip(begin_list, end_list):
        yield begin, end


def make_pos(begin, end):
    return np.arange(begin, end) + 1


def check_outfile(outfile):
    if pathlib.Path(outfile).exists():
        logger.error(
            f"{outfile} exists! Make sure to provide "
            "a non-existing output file name"
        )
        sys.exit()


def to_str(array):
    """Convert array of arrays to string"""
    return "\n".join("\t".join(list(map(str, x))) for x in array)


def parse_region(region):
    m = re.match(r"^(?P<chrom>[\w]+):?(?P<begin>\d+)?-?(?P<end>\d+)?", region)
    if m is None:
        logger.error(f"Malformatted regions string: {region}")
        sys.exit(1)

    return m.groups()


class D4Iterator:
    """Iterate over multiple d4 paths"""

    def __init__(self, path, chunk_size=10000, regions=None, concat=False):
        self._fh = [pyd4.D4File(x) for x in tqdm(path)]
        self._index = len(self._fh)
        self._chunk_size = chunk_size
        if regions is None:
            if concat:
                self._chroms = [x for fh in self._fh for x in fh.chroms()]
            else:
                self._chroms = self._fh[0].chroms()
        else:
            self._chroms = []
            for chrom_name, end in self._fh[0].chroms():
                reg = [str(chrom_name), 0, end]
                if regions.T.isin(reg).all().any():
                    self._chroms.append((chrom_name, end))

    @property
    def chroms(self):
        return self._chroms

    @property
    def chunk_size(self):
        return self._chunk_size

    @property
    def writer(self):
        return self._writer

    @writer.setter
    def writer(self, fn):
        self._writer = (
            pyd4.D4Builder(str(fn)).add_chroms(self.chroms).get_writer()
        )

    def __iter__(self):
        return self

    def __next__(self):
        if self._index == 0:
            self._index = len(self._fh)
            raise StopIteration
        self._index = self._index - 1
        return self._fh[self._index]

    def iter_chroms(self):
        for chrom_name, end in (pbar := tqdm(self.chroms)):
            pbar.set_description(f"processing chromosome {chrom_name}")
            yield chrom_name, 0, end

    def iter_chunks(self, chrom_name, begin, end):
        for rbegin, rend in (
            pbar := tqdm(make_chunks(begin, end, self.chunk_size))
        ):
            rname = f"{chrom_name}:{rbegin}-{rend}"
            pbar.set_description(f"processing region {rname}")
            yield rname

    def process_region_chunk(self, rname):
        for i, track in (pbar := tqdm(enumerate(self))):
            pbar.set_description(f"processing track {i}")
            yield i, track.load_to_np(rname)

    def sum(self, chrom_name, begin, end):  # noqa: A003
        """Sum tracks over a chromosome region"""

        def _sum_region_chunk(rname):
            for i, y in self.process_region_chunk(rname):
                if i == 0:
                    x = y
                else:
                    x = x + y
            return x

        for j, rname in enumerate(self.iter_chunks(chrom_name, begin, end)):
            if j == 0:
                y = _sum_region_chunk(rname)
            else:
                y = np.append(y, _sum_region_chunk(rname))
        return y

    def count(self, chrom_name, begin, end, *, lower=0, upper=np.inf):
        """Count tracks whose values lie in a given range over a given
        region"""

        def _count_region_chunk(rname):
            for i, y in self.process_region_chunk(rname):
                if i == 0:
                    x = ((y > lower) & (y < upper)).astype(int)
                else:
                    x = x + ((y > lower) & (y < upper)).astype(int)
            return x

        for j, rname in enumerate(self.iter_chunks(chrom_name, begin, end)):
            if j == 0:
                y = _count_region_chunk(rname)
            else:
                y = np.append(y, _count_region_chunk(rname))
        return y


# FIXME: regions should be parameter to *all* functions
@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__


@cli.command()
@click.argument("path", nargs=-1, type=click.Path(exists=True))
@click.argument("outfile", type=click.Path(exists=False))
@click.option(
    "--chunk-size", help="region chunk size", default=1000000, type=int
)
@click.option("--regions", "-R", help="region bed file")
@click.pass_context
def sum(  # noqa: A001
    ctx,
    path,
    outfile,
    chunk_size,
    regions,
):
    """Sum d4 files.

    Sum first track from multiple d4 files to a single-track file.

    """
    logger.info(
        f"Running {ctx.find_root().info_name} version " + ctx.obj["VERSION"]
    )
    check_outfile(outfile)

    bed = None
    if regions is not None:
        bed = pd.read_table(
            regions,
            names=["chrom", "begin", "end"],
            usecols=[0, 1, 2],
            header=None,
        )

    d4fh = D4Iterator(path, chunk_size=chunk_size, regions=bed)
    d4fh.writer = outfile
    for chrom_name, begin, end in d4fh.iter_chroms():
        y = d4fh.sum(chrom_name, begin, end)
        d4fh.writer.write_np_array(chrom_name, 0, y)
    d4fh.writer.close()


@cli.command()
@click.argument("path", nargs=-1, type=click.Path(exists=True))
@click.argument("outfile", type=click.Path(exists=False))
@click.option("--chunk-size", help="region chunk size", default=1000000)
@click.option("--min-coverage", help="minimum coverage", default=0, type=int)
@click.option("--max-coverage", help="maximum coverage", type=int)
@click.option("--regions", "-R", help="region bed file")
@click.pass_context
def count(ctx, path, outfile, chunk_size, min_coverage, max_coverage, regions):
    """Count coverages in files.

    Count the number of tracks that have coverages in a given range.
    """
    logger.info(
        f"Running {ctx.find_root().info_name} version " + ctx.obj["VERSION"]
    )
    check_outfile(outfile)

    bed = None
    if regions is not None:
        bed = pd.read_table(
            regions,
            names=["chrom", "begin", "end"],
            usecols=[0, 1, 2],
            header=None,
        )

    if max_coverage is None:
        max_coverage = np.inf
    d4fh = D4Iterator(path, chunk_size=chunk_size, regions=bed)
    d4fh.writer = outfile
    for chrom_name, begin, end in d4fh.iter_chroms():
        y = d4fh.count(
            chrom_name, begin, end, lower=min_coverage, upper=max_coverage
        )
        d4fh.writer.write_np_array(chrom_name, 0, y)
    d4fh.writer.close()


@cli.command()
@click.argument("path", type=click.Path(exists=True))
@click.argument("regions", type=str, nargs=-1, default=None, required=False)
@click.option("--outfile", "-o", default=sys.stdout, type=click.File("w"))
@click.option("--header", "-h", help="Include header", is_flag=True)
@click.option(
    "--window-size", "-w", help="window size", type=int, default=None
)
@click.option("--step-size", "-s", help="step size", type=int, default=None)
@click.option("--exclude", "-e", help="regexp to exclude regions", type=str)
@click.pass_context
def view(ctx, path, regions, outfile, header, window_size, step_size, exclude):
    """View d4 file enumerating all positions.

    View a d4 file enumerating all positions or by windows.

    """
    fh = pyd4.D4File(str(path))
    chrom_sizes = dict(fh.chroms())
    if header:
        outfile.write("#CHROM\tPOS\tVALUE\n")
    if len(regions) == 0:
        regions = fh.chrom_names()
    for reg in regions:
        logger.debug(f"Processing region {reg}")
        chrom_name, begin, end = parse_region(reg)
        try:
            y = fh.load_to_np(reg)
        except KeyError:
            logger.warn(f"No such region {reg}")
            continue
        if exclude is not None:
            if re.search(rf"{exclude}", reg):
                logger.debug(f"Skipping region {reg}")
                continue
        if window_size is not None:
            if window_size > len(y):
                window_size = len(y)
            if begin is None:
                begin = 0
            else:
                begin = int(begin)
            if end is None:
                end = chrom_sizes[chrom_name]
            else:
                end = int(end)
            window_starts = np.arange(begin, end, window_size)
            window_stops = np.arange(
                begin + window_size, end + window_size, window_size
            )
            window_stops[-1] = min(window_stops[-1], end)
            if step_size is None:
                step_size = window_size
            windows = sliding_window_view(y, window_size)
            data = windows[::step_size, :]
            imax = window_size * len(data)
            stat = data.mean(axis=1)
            if imax < len(y):
                stat = np.append(stat, y[imax:].mean())
            for wstart, wstop, val in zip(window_starts, window_stops, stat):
                outfile.write(f"{chrom_name}\t{wstart}\t{wstop}\t{val}\n")
        else:
            for _, pos, val in fh.enumerate_values(
                chrom_name, 0, chrom_sizes[chrom_name]
            ):
                vout = "\t".join([str(int(x)) for x in val])
                outfile.write(vout)


@cli.command()
@click.argument("path", type=click.Path(exists=True))
@click.argument("outfile", type=click.Path(exists=False))
@click.option("--chunk-size", help="region chunk size", default=1000000)
@click.option("--lower", help="lower bound", default=0, type=int)
@click.option("--upper", help="upper bound", default=np.inf, type=int)
@click.option("--regions", "-R", help="region bed file")
@click.pass_context
def filter(ctx, path, outfile, chunk_size, lower, upper, regions):  # noqa: A001
    """Create boolean output from filter.

    Filter file on value range to generate boolean-encoded output
    (i.e. a sequence mask).

    """
    logger.info(
        f"Running {ctx.find_root().info_name} version " + ctx.obj["VERSION"]
    )
    check_outfile(outfile)

    bed = None
    if regions is not None:
        bed = pd.read_table(
            regions,
            names=["chrom", "begin", "end"],
            usecols=[0, 1, 2],
            header=None,
        )

    fh = pyd4.D4File(str(path))
    for region in tqdm(bed.itertuples()):
        x = fh[f"{region.chrom}:{region.begin}-{region.end}"]
        print(x[0:10])
    # d4fh = D4Iterator([path], chunk_size=chunk_size, regions=bed)
    # #output_writer =
    # #pyd4.D4Builder(str(outfile)).add_chroms(fh.chroms()).get_writer()
    # d4fh.writer = outfile
    # for chrom_name, begin, end in d4fh.iter_chroms():
    #     y = (d4fh.count(chrom_name, begin, end, lower=lower,
    #     upper=upper) > 0).astype(int)
    #     d4fh.writer.write_np_array(chrom_name, 0, y)
    # d4fh.writer.close()
    # # for chrom_name in tqdm(fh.chrom_names()):
    # #     x = ((fh[chrom_name] > lower) & (fh[chrom_name] <
    # #     upper)).astype(int)
    # #     output_writer.write_np_array(chrom_name, 0, x)
    # # output_writer.close()


@cli.command()
@click.argument("path", nargs=-1, type=click.Path(exists=True))
@click.argument("outfile", type=click.Path(exists=False))
@click.pass_context
def concat(ctx, path, outfile):
    """Concatenate d4 output.

    Concatenate d4 files given by PATH to OUTFILE.
    """
    logger.info(
        f"Running {ctx.find_root().info_name} version " + ctx.obj["VERSION"]
    )
    check_outfile(outfile)

    d4fh = D4Iterator(path, concat=True)
    d4fh.writer = outfile
    for fn, fh in (pbar := tqdm(zip(path, d4fh))):
        pbar.set_description(f"processing file {fn}")
        for chrom_name in (pinnerbar := tqdm(fh.chrom_names())):
            pinnerbar.set_description(f"processing chromosome {chrom_name}")
            y = fh[chrom_name]
            d4fh.writer.write_np_array(chrom_name, 0, y)
    d4fh.writer.close()


# Best? Run d4tools on region, bedtools to summarize over feature.
# Need feature definition files as input.
# @cli.command()
# @click.argument("path", nargs=-1, type=click.Path(exists=True))
# @click.argument("outfile", type=click.Path(exists=False))
# @click.argument("regions", type=click.Path(exists=True))
# @click.pass_context
# def matrix(ctx, path, outfile):
#     """Create matrix from d4 files.

#     Create a matrix from d4 files given by PATH to OUTFILE. Each row
#     begins with a BED region definition followed by the values for
#     each track. The output is plain text as d4 format does not
#     automatically reduce file size."""
#     logger.info(
#         f"Running {ctx.find_root().info_name} version " + ctx.obj["VERSION"]
#     )
#     check_outfile(outfile)

#     d4fh = D4Iterator(path)


# @cli.command()
# @click.argument("path", type=click.Path(exists=True))
# @click.argument("outfile", type=click.Path(exists=False))
# @click.option(
#     "--num-sites", "-n", help="number of sites", type=int, default=10_000
# )
# @click.option("--regions", "-R", help="region bed file")
# @click.option("--chunk-size", help="region chunk size", default=1000000)
# @click.pass_context
# def sample(ctx, path, outfile, num_sites, regions, chunk_size):
#     """Sample d4 file.

#     Sample d4 file to a given number of sites, possibly restricted to
#     pre-defined regions.

#     """

#     logger.info(
#         f"Running {ctx.find_root().info_name} version " + ctx.obj["VERSION"]
#     )
#     check_outfile(outfile)

#     bed = None
#     if regions is not None:
#         bed = pd.read_table(
#             regions,
#             names=["chrom", "begin", "end"],
#             usecols=[0, 1, 2],
#             header=None,
#         )

#     if max_coverage is None:
#         max_coverage = np.inf
#     d4fh = D4Iterator(path, chunk_size=chunk_size, regions=bed)

#     fh = pyd4.D4File(str(path))
#     chrom_sizes = dict(fh.chroms())
#     for chrom_name, end in chrom_sizes.items():
#         y = fh.load_to_np(f"{chrom_name}:0-{end}")
#         if num_sites < len(y):
#             idx = np.random.choice(len(y), num_sites, replace=False)
#             y = y[idx]
#         fh.write_np_array(chrom_name, 0, y)
#     fh.close()
