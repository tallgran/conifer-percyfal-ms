"""Summarize diversity statistics

Summarize diversity statistics in INFILE and write to OUTFILE. The
INFILE is a file with position and score information, and could be a
tab-separated file as produced by vcftools, or a BEDGRAPH file. The
output format is compressed BEDGRAPH.

In addition to the INFILE, a GENOME file in FASTA INDEX format is
required that provides information on chromosome names and sizes.

Genomic sites can be masked by ACCESSIBILITY MASK files, which are BED
files that define accessible regions. The regions are converted into
0/1-numpy arrays. More than one mask can be supplied, in which case
the masks will be combined. The number of accessible sites will be
summarized over the genome and/or regions that are being analyzed.

In addition, FEATURE definition files are BED-formatted files that can
restrict analyses to regions or features of interest.

Finally, it is possible to define a WINDOW size or a REGION file that
determines regions over which statistics will be summarized. Note that
this differs from the FEATURE definition file in that a feature may be
split between WINDOWs. The REGION file should be supplied in BED
format.

"""

import concurrent.futures
import logging
import sys
from threading import BoundedSemaphore

import click
import numpy as np
import pandas as pd
from pympler import asizeof
from tqdm import tqdm

from conifer.io.bed import read_fai_genomefile
from conifer.io.csv import BufferedCsvReader
from conifer.io.files import size
from conifer.io.vcftools import BufferedVcftoolsReader
from conifer.model.region import (
    AnalysisChunk,
    Bed4Region,
    Region,
    Score,
    ScoreVector,
    calculate_summary_stats,
    dataframe_to_region_list,
    subset_dataframe_by_region,
)

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)


np.seterr(divide="ignore", invalid="ignore")


def compress_data_chunk(df, log_level=3):
    """Compress data chunk"""
    if log_level >= 3:
        s1, s1hum = size(df, False), size(df)
    df["CHROM"] = pd.Categorical(
        df["CHROM"], categories=df["CHROM"].unique(), ordered=True
    )
    if log_level >= 3:
        s2, s2hum = size(df, False), size(df)
    logger.debug(
        "Compressed data chunk: %s -> %s (%.2f%% reduction)",
        s1hum,
        s2hum,
        100 * (1 - s2 / s1),
    )


def handle_termination_signal(sig, frame):  # pylint: disable=unused-argument
    """Handle termination signal"""
    sys.exit(0)


def read_genomefile(file) -> pd.DataFrame:
    """Read fasta index genome file"""
    df = read_fai_genomefile(file)
    df = df[["chrom", "begin", "end"]]
    if df is not None:
        data = {}
        for group, group_data in df.groupby("chrom"):
            name = str(Region(*group_data.iloc[0]))
            data[group] = Bed4Region(*group_data.iloc[0], name=name)
        return data
    return None


def read_bedfile(file, filetype="bed") -> pd.DataFrame:
    """Load a single BED file"""
    logger.info("Processing BED file %s", file.name)
    kwargs = {
        "names": ["chrom", "begin", "end"],
        "dtype": {"chrom": "category", "begin": np.int64, "end": np.int64},
    }
    if filetype == "bed4":
        kwargs["names"].append("name")
        kwargs["dtype"]["name"] = str
    reader = BufferedCsvReader(file.name, chunksize=None, **kwargs)
    return next(reader)


def load_bedfiles(infiles, labels=None) -> dict[pd.DataFrame]:
    """Load BED files. Return dict with filename to data mapping."""
    logger.info("Loading bedfiles %s", ",".join([fn.name for fn in infiles]))
    results = {}
    try:
        labels = labels.split(",")
    except AttributeError:
        labels = [fn.name for fn in infiles]
    for fn, name in zip(infiles, labels):
        data = read_bedfile(fn)
        results[name] = data
    return results


def make_windows(chrom, window_size, step_size=None) -> list[Region]:
    """Make windows sized window_size over a sequence of length
    sequence_length"""
    if step_size is None:
        step_size = window_size
    begin = np.arange(0, len(chrom), window_size)
    end = np.append(begin[1:], len(chrom))
    result = [
        Bed4Region(chrom.chrom, coords[0], coords[1], f"w{i}")
        for i, coords in enumerate(zip(begin, end))
    ]
    return result


def make_regions(
    chromosomes, regionfile=None, window_size=None
) -> pd.DataFrame:
    """Define regions over which to summarize results. Defaults to
    entire chromosomes."""
    logger.info("Defining regions")
    if (window_size is None) and (regionfile is None):
        return {"genome": pd.DataFrame(chromosomes.values())}
    if window_size is not None:
        windows = []
        for chrom in chromosomes.values():
            windows.extend(make_windows(chrom, window_size))
        return {f"w{window_size}": pd.DataFrame(windows)}
    if regionfile is not None:
        custom = read_bedfile(regionfile, filetype="bed4")
        return {"custom": custom}
    return None


def process_chunk(args) -> pd.DataFrame:  # pylint: disable=too-many-locals
    """Process analysis chunk"""
    _, chunk = args
    logger.info("Processing chunk %s", chunk)
    results = []
    for rname, regions in chunk.region.items():
        logger.debug("Region name: %s", rname)
        for reg in tqdm(regions, desc=f"Regions in {chunk}"):
            logger.info(
                "Subsetting chunk %s to %i - %i",
                str(chunk),
                reg.begin,
                reg.end,
            )
            item = slice(reg.begin - chunk.begin, reg.end - chunk.begin, None)
            rchunk = chunk[item]

            score_vector = ScoreVector(
                chrom=rchunk.chrom,
                begin=rchunk.begin,
                end=rchunk.end,
                name=rchunk.name,
                score_vector=dataframe_to_region_list(
                    rchunk.score,
                    cls=Score,
                    strand=".",
                    name=rname,
                ),
            )

            for ftname in rchunk.feature:
                (
                    score,
                    total_score,
                    site_score,
                    n_sites,
                    n_accessible,
                    n_segregating_sites,
                    n_segregating_sites_accessible,
                    theta,
                    tajima,
                ) = calculate_summary_stats(
                    reg,
                    rchunk.amask,
                    rchunk.feature[ftname],
                    score_vector.score,
                    score_vector.n_chr,
                    score_vector.score_pos,
                )
                results.append(
                    [
                        reg.chrom,
                        reg.begin,
                        reg.end,
                        reg.name,
                        ftname,
                        chunk.name,
                        score,
                        total_score,
                        site_score,
                        n_sites,
                        n_accessible,
                        n_segregating_sites,
                        n_segregating_sites_accessible,
                        theta,
                        tajima,
                    ]
                )
    return pd.DataFrame(results)


class MaxQueuePool:
    """This Class wraps a concurrent.futures.Executor limiting the
    size of its task queue.

    If `max_queue_size` tasks are submitted, the next call to submit
    will block until a previously submitted one is completed.

    cf https://gist.github.com/noxdafox/4150eff0059ea43f6adbdd66e5d5e87e
    """

    def __init__(self, executor, *, max_queue_size, max_workers=None):
        logger.info(
            "Initializing queue with %i queue slots, %i workers",
            max_queue_size,
            max_workers,
        )
        self.pool = executor(max_workers=max_workers)
        self.pool_queue = BoundedSemaphore(max_queue_size)

    def submit(self, fn, *args, **kwargs):
        """Submit a new task to the pool. This will block if the queue
        is full"""
        self.pool_queue.acquire()  # pylint: disable=consider-using-with
        future = self.pool.submit(fn, *args, **kwargs)
        future.add_done_callback(self.pool_queue_callback)

        return future

    def pool_queue_callback(self, _):
        """Called when a future is done. Releases one queue slot."""
        self.pool_queue.release()


@click.command(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.argument("infile", type=click.File("r"))
@click.argument("countsfile", type=click.File("r"))
@click.argument("genomefile", type=click.File("r"))
@click.option("--outfile", "-o", type=click.File("w"), default=sys.stdout)
@click.option(
    "--accessibility-mask",
    "-a",
    "maskfiles",
    help="Accessibility mask file in BED format",
    multiple=True,
    type=click.File("r"),
)
@click.option(
    "--feature",
    "-f",
    "featurefiles",
    help="Feature definition files in BED format",
    multiple=True,
    type=click.File("r"),
)
@click.option("--window-size", "-w", help="Window size", type=int)
@click.option(
    "--region",
    "-r",
    "regionfile",
    help="Region definition file in BED format",
    type=click.File("r"),
)
@click.option(
    "--mask-labels",
    help="Comma-separated string of accessibility mask labels",
    type=str,
)
@click.option(
    "--feature-labels",
    help="Comma-separated string of feature labels",
    type=str,
)
@click.option(
    "--chunksize",
    help="Chunksize (lines) to read from infile",
    type=int,
    default=1000000,
)
@click.option(
    "--region-chunksize",
    help="Size (bp) of processed region chunks.",
    type=int,
    default=1000000,
)
@click.option(
    "--cores",
    "-j",
    help="Number of CPUs",
    type=int,
    default=1,
)
@click.option(
    "--statistic",
    "-s",
    help="Label for statistic in output",
    type=str,
    default="pi",
)
@click.option(
    "--verbose",
    "-v",
    help="Set the verbosity level",
    count=True,
)
@click.pass_context
def cli(
    ctx,
    infile,
    countsfile,
    genomefile,
    maskfiles,
    featurefiles,
    outfile,
    window_size,
    regionfile,
    mask_labels,
    feature_labels,
    chunksize,
    region_chunksize,
    cores,
    statistic,
    verbose,
):  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    """Summarize_diversity cli"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.info("Running %s version %s", __name__, __version__)

    log_level = max(3 - verbose, 0) * 10
    logging.basicConfig(
        level=log_level,
        format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
    )

    if (window_size is not None) and (regionfile is not None):
        raise click.UsageError(
            "--window-size and --region cannot be supplied simultaneously"
        )

    chromosomes = read_genomefile(genomefile)

    mask = load_bedfiles(maskfiles, mask_labels)

    features = load_bedfiles(featurefiles, feature_labels)

    if not features:
        features = make_regions(chromosomes, None, None)

    regions = make_regions(chromosomes, regionfile, window_size)
    reader = BufferedVcftoolsReader(
        infile.name, chunksize=chunksize, position_chunksize=region_chunksize
    )
    counts_reader = BufferedVcftoolsReader(
        countsfile.name,
        chunksize=chunksize,
        sep="\t",
        position_chunksize=region_chunksize,
        names=["CHROM", "POS", "N_CHR"],
        usecols=[0, 1, 3],
    )

    def _generate_analysis_chunks():
        old_chrom = None
        begin = 0
        for data, count_data in zip(
            (pbar := tqdm(reader, disable=False)), counts_reader
        ):
            chrom, df = data
            _, counts = count_data
            df = pd.merge(df, counts)
            compress_data_chunk(df, log_level)
            df.columns = ["chrom", "end", "score", "n_chr"]
            df["begin"] = df["end"] - 1
            df = df[["chrom", "begin", "end", "score", "n_chr"]]
            if chrom != old_chrom:
                old_chrom = chrom
                begin = 0
            else:
                begin += region_chunksize

            end = min(
                begin + region_chunksize,
                chromosomes[chrom].end,
            )

            chunk = AnalysisChunk(
                chrom=chrom, begin=begin, end=end, name=statistic
            )
            for label, m in mask.items():
                submask = subset_dataframe_by_region(m, chunk)
                chunk.add_mask(label, submask, invert=True)

            for label, ft in features.items():
                subfeature = subset_dataframe_by_region(ft, chunk)
                chunk.add_feature(label, subfeature, invert=True)

            for label, reg in regions.items():
                subregion = subset_dataframe_by_region(reg, chunk)
                chunk.add_region(label, subregion)
            chunk.score = df
            pbar.set_description(
                (
                    f"Generated {chunk} chunk "
                    f"({asizeof.asizeof(chunk) / 1024.0 / 1000.0:.2f}M, "
                    f"{len(chunk)}bp)"
                ),
                refresh=False,
            )
            if (len(chunk) == 0) or (len(chunk) > region_chunksize):
                raise ValueError(chunk.score)
            yield chrom, chunk

    generator = _generate_analysis_chunks()
    futures = []
    pool = MaxQueuePool(
        concurrent.futures.ProcessPoolExecutor,
        max_workers=cores,
        max_queue_size=int(2 * cores),
    )
    for args in generator:
        futures.append(pool.submit(process_chunk, args))

    result = pd.concat([x.result() for x in futures])
    result.columns = [
        "chrom",
        "begin",
        "end",
        "region",
        "feature",
        "statistic",
        "score",
        "total_score",
        "site_score",
        "n_sites",
        "n_accessible",
        "n_segregating_sites",
        "n_segregating_sites_accessible",
        "theta",
        "tajimasD",
    ]

    result["chrom"] = pd.Categorical(
        result["chrom"], categories=result["chrom"].unique(), ordered=True
    )
    result = result.sort_values(by=["chrom", "begin", "region", "feature"])
    result.to_csv(outfile, index=False)
    outfile.close()
