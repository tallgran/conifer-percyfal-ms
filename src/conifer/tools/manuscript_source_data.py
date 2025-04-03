"""manuscript_source_data

manuscript_source_data documentation.
"""

import csv
import logging
import os
import re
import sys
from collections import namedtuple
from pathlib import Path

import click
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from scipy.stats import bootstrap
from tqdm import tqdm

from conifer.options import verbose_option

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)

REFERENCE_SAMPLES = {
    "diploid": "Diploid",
    "haploid_ERX242653": "Haploid_ERX242654",
    "Engelmannii_SRX1530215": "P_engelmannii",
    "Glauca_SRX160982": "P_glauca",
    "Sitchensis_SRX2354260": "P_sitchensis",
}

HARMONIC_NUMBER_MAX = 2500

HARMONIC_NUMBER = np.zeros(HARMONIC_NUMBER_MAX + 1)
HARMONIC_NUMBER_SQUARED = np.zeros(HARMONIC_NUMBER_MAX + 1)
HARMONIC_NUMBER[0] = np.nan
HARMONIC_NUMBER_SQUARED[0] = np.nan
for ii in np.arange(1, HARMONIC_NUMBER_MAX):
    HARMONIC_NUMBER[ii] = np.sum(1 / np.arange(1, ii))
    HARMONIC_NUMBER_SQUARED[ii] = np.sum(1 / (np.arange(1, ii) ** 2))

EULER_MASCHERONI = 0.57721566490153286060651209008240243104215933593992

# The call depth is estimated as mean from the count masks
THETA_CALL_DEPTH = {
    "P.abies": 1529,
    "highcov": 567,
    "north": 50,
    "south": 49,
}

SAMPLESETS = [
    "P.abies",
    "highcov",
    "north",
    "south",
    # "recombination",
]
SAMPLESET_DESC = pd.DataFrame(
    {
        "Sampleset": SAMPLESETS,
        "n": [0] * len(SAMPLESETS),
        "Description": [
            "All samples",
            "High coverage samples",
            "Northern samples",
            "Southern samples",
            # "High coverage samples, North",
            # "High coverage samples, South",
            # "Samples used for recombination analyses",
            # "Rangewide samples",
        ],
    }
).set_index("Sampleset")

FEATURES = [
    "CDS",
    "UTR",
    "intron",
    "TE.gene",
    "pseudogene",
    "intergenic",
    "genome",
]
DataPath = namedtuple("DataPath", ["source_data", "figshare"])


def get_feature_size(ds, features):
    """Return the size of a feature"""
    results = []
    for g, d in ds.groupby("feature"):
        n_sites = d["n_sites"].sum()
        size = convert_to_si_suffix(n_sites)
        results.append((g, n_sites, size))
    data = pd.DataFrame(results, columns=["feature", "n_sites", "size"])
    data["feature"] = pd.Categorical(
        data["feature"], categories=features, ordered=True
    )
    data = data.sort_values("feature")
    return data


def convert_to_si_suffix(x):
    """Convert number to SI suffix"""
    suffixes = ["", "k", "M", "G"]
    order_magnitude = np.floor(len(str(abs(x))) / 3)
    sfx = suffixes[int(order_magnitude)]
    scaled_value = x / 10 ** (3 * order_magnitude)
    return f"{scaled_value:.1f}{sfx}bp"


def data_path_callback(ctx, param, value):  # pylint: disable=unused-argument
    """Set output file path for source data and fishare data"""
    if value is None:
        value = Path(os.environ.get("MANUSCRIPT_ROOT"))
    else:
        value = Path(value)
    path = DataPath(
        value / "mainGenome2024" / "source_data",
        value / "mainGenome2024" / "figshare",
    )
    if not path.source_data.exists():
        logger.info("Creating source data directory: %s", path.source_data)
        path.source_data.mkdir(parents=True)
    if not path.figshare.exists():
        logger.info("Creating figshare data directory: %s", path.figshare)
        path.figshare.mkdir(parents=True)
    logger.info("Data path: %s", path)
    return path


data_path_opt = click.option(
    "--data-path",
    type=str,
    help="Set the data path",
    default=None,
    callback=data_path_callback,
    expose_value=True,
)


def paccessible_opt(default=0.5):
    return click.option(
        "--paccessible",
        type=float,
        default=default,
        help="Minimum proportion of accessible sites",
    )


def total_sites_opt(default=1000):
    return click.option(
        "--total-sites",
        type=int,
        default=default,
        help="Minimum number of total sites in a window",
    )


def samplesets_opt(default):
    """Samplesets option"""

    def samplesets_callback(ctx, param, value):  # pylint: disable=unused-argument
        """Samplesets callback"""
        logger.info("Sample sets: %s", value)
        return value

    multiple = isinstance(default, list)
    return click.option(
        "--samplesets",
        type=click.Choice(SAMPLESETS),
        multiple=multiple,
        default=default,
        callback=samplesets_callback,
        help="Sample sets to include",
    )


def features_opt(default):
    """Features option"""

    if default is None:
        default = FEATURES

    def features_callback(ctx, param, value):  # pylint: disable=unused-argument
        """Features callback"""
        logger.info("Features: %s", value)
        return value

    return click.option(
        "--features",
        type=click.Choice(default),
        multiple=True,
        default=default,
        callback=features_callback,
        help="Features to include",
    )


def regex_opt(default):
    """Regex option"""
    return click.option(
        "--regex", type=str, help="Input file regex", default=default
    )


def output_file_regex_opt(default):
    """Output file regex option"""
    return click.option(
        "--output-file-regex",
        type=str,
        help="Output file regex",
        default=default,
    )


def window_size_opt(default):
    """Window size option"""

    def window_size_callback(ctx, param, value):  # pylint: disable=unused-argument
        """Window size callback"""
        logger.info("Window size: %s", value)
        return int(value)

    multiple = isinstance(default, list)
    return click.option(
        "--window-size",
        type=click.Choice(["1000", "50000", "100000", "1000000", "10000000"]),
        multiple=multiple,
        callback=window_size_callback,
        default=default,
    )


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """manuscript_source_data docstring"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.info("Running %s version %s", __name__, __version__)


@cli.command()
@data_path_opt
@samplesets_opt(default=SAMPLESETS)
@regex_opt("resources/samplesets/samples-{sampleset}.tsv")
@verbose_option()
def make_samplesets_table(data_path, samplesets, regex):
    """Make samplesets table"""
    logger.info("Making samplesets table")
    infiles = [regex.format(sampleset=ss) for ss in samplesets]
    for ss, infile in zip(samplesets, infiles):
        data = pd.read_csv(infile, sep="\t")
        SAMPLESET_DESC.loc[ss, "n"] = data.shape[0]
    output_file = data_path.source_data / "tbl-samplesets.csv"
    logger.info("Writing %s", output_file)
    SAMPLESET_DESC.to_csv(
        output_file, sep=",", index=True, quoting=csv.QUOTE_NONNUMERIC
    )


@cli.command()
@data_path_opt
@click.option(
    "--samplesheet",
    type=click.Path(exists=True),
    help="Samplesheet",
    default="resources/samplesheet.tsv",
)
@click.option(
    "--sample-xml",
    type=click.Path(exists=True),
    help="Submission XML",
    default="data/resources/upsc/UPSC-0232.Sample.xml",
)
@click.option(
    "--receipt-xml",
    type=click.Path(exists=True),
    help="Submission XML receipt",
    default="data/resources/upsc/UPSC-0232.Receipt.xml",
)
@click.option(
    "--samtools-stats-regex",
    type=str,
    help="Regex for samtools stats files",
    default="results/samtools_stats/{sample}.txt",
)
@click.option(
    "--reference-fai",
    type=click.Path(exists=True),
    help="Reference fasta index",
    default="data/resources/pabies-2.0.fa.fai",
)
@click.option(
    "--show",
    is_flag=True,
    help="Show data",
)
@click.option(
    "--full",
    is_flag=True,
    help="Full output",
)
@verbose_option()
def make_samplesheet_table(
    data_path,
    samplesheet,
    sample_xml,
    receipt_xml,
    samtools_stats_regex,
    reference_fai,
    show,
    full,
):
    """Make samplesheet table for manuscript"""
    ss = pd.read_csv(samplesheet, sep="\t")
    ss["VCF-ID"] = ss["SM"].copy()
    ss["SM"] = ss["SM"].str.replace(r"-", "_", regex=True)

    with open(sample_xml, "r") as f:
        data = f.read()
    soup = BeautifulSoup(data, "xml")
    title2alias = []
    for sample in soup.find_all("SAMPLE"):
        title2alias.append((sample.find("TITLE").text, sample["alias"]))

    df = pd.DataFrame(title2alias, columns=["title", "Alias"])
    df["ena-sample-accession"] = None
    df["ena-biosample"] = None

    receipt = BeautifulSoup(open(receipt_xml, "r"), "xml")

    for sample in receipt.find_all("SAMPLE"):
        df.loc[df["Alias"] == sample["alias"], "ena-sample-accession"] = (
            sample["accession"]
        )
        for ext_id in sample.find_all("EXT_ID"):
            if ext_id["type"] == "biosample":
                df.loc[df["Alias"] == sample["alias"], "ena-biosample"] = (
                    ext_id["accession"]
                )

    dfout = ss.merge(df, left_on="SM", right_on="title", how="left")
    dfout["Biosample"] = dfout["Biosample"].replace("", np.nan)
    dfout["Biosample"] = dfout["Biosample"].fillna(dfout["ena-biosample"])
    dfout["Experiment-ID"] = dfout["Experiment-ID"].replace("", np.nan)
    dfout["ENA-project-ID"] = dfout["ENA-project-ID"].fillna("PRJEB85275")

    # Remove ONT samples
    ix = dfout["Note"].str.contains("ONT", na=False)
    dfout = dfout[~ix]

    dfout["mapped_dedup_sequences"] = None
    dfout["bases_mapped"] = None
    for i, row in tqdm(dfout.iterrows()):
        fn = samtools_stats_regex.format(sample=row["SM"])
        if not os.path.exists(fn):
            fn = samtools_stats_regex.format(sample=row["Scilife-ID"])
        logger.info("Reading %s", fn)
        try:
            with open(fn, "r") as f:
                data = "".join(f.read())
                sequences = re.findall(r"SN\s+sequences:\s+(\d+)", data)
                bases_mapped = re.findall(
                    r"SN\s+bases mapped \(cigar\):\s+(\d+)", data
                )
            dfout.loc[i, "mapped_dedup_sequences"] = sequences[0]
            dfout.loc[i, "bases_mapped"] = bases_mapped[0]
        except FileNotFoundError:
            logger.warning(
                "File not found: %s",
                samtools_stats_regex.format(sample=row["SM"]),
            )
        except IndexError:
            logger.warning("No sequences found for %s", row["SM"])

    fai = pd.read_csv(reference_fai, sep="\t", header=None)
    fai.columns = ["chrom", "length", "offset", "linebases", "linewidth"]
    genome_length = fai["length"].sum()
    dfout["mapping_coverage"] = (
        dfout["bases_mapped"].astype(float) / genome_length
    )
    dfout["Scilife-ID"] = ss["SM"].copy()
    dfout["Scilife-ID"] = dfout["Scilife-ID"].replace(REFERENCE_SAMPLES)
    dfout["VCF-ID"] = dfout["VCF-ID"].replace(REFERENCE_SAMPLES)
    output_file = data_path.source_data / "tbl-samplesheet-full.csv"
    if not full:
        output_file = data_path.source_data / "tbl-samplesheet.csv"
        outcols = [
            "Biosample",
            "Lab-ID",
            "VCF-ID",
            "PopulationType",
            "Latitude",
            "Longitude",
            "mapping_coverage",
        ]
        outcol_names = [
            "ENA biosample accession",
            "Internal ID",
            "VCF ID",
            "Population type",
            "Latituted",
            "Longitude",
            "Mapping coverage",
        ]
        dfout = dfout[outcols]
        dfout.columns = outcol_names
        # dfout = dfout[dfout["ENA biosample accession"].isin(df["biosample"])]
    if show:
        dfout.to_csv(
            sys.stdout, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
        )
    else:
        logger.info("Writing to %s", output_file)
        dfout.to_csv(
            output_file, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
        )


@cli.command()
@data_path_opt
@samplesets_opt(default=SAMPLESETS[0:6])
@features_opt(FEATURES)
@output_file_regex_opt("fig-individual-coverage-{sampleset}.csv")
@regex_opt(
    "data/aggregate_coverage/{sampleset}/MQ10/count_ge3.p100.d4.{feature}.hist"
)
@click.option(
    "--sample-size", type=int, default=10000, help="Number of samples to draw"
)
@verbose_option()
def compile_coverage(
    data_path, output_file_regex, samplesets, features, regex, sample_size
):
    """Compile individual coverage data.

    Compile output of mosdepth histogram files located in
    data/aggregate_coverage. The histograms are used to generate a
    weighted sample with replacement.

    """
    logger.info("Compiling coverage data")

    for ss in samplesets:
        for ft in features:
            infile = regex.format(sampleset=ss, feature=ft)
            try:
                df = pd.read_csv(
                    infile,
                    sep="\t",
                    header=None,
                    names=["x", "count"],
                    dtype={"x": str, "count": np.int32},
                )
            except FileNotFoundError:
                logger.error("File not found: %s", infile)
                continue
            df = df.replace("<0", -1)
            df = df.replace(r">(\d+)", r"\1", regex=True)
            df["x"] = df["x"].astype(np.int32)
            df["sampleset"] = ss
            df["feature"] = ft
            if "data" not in locals():
                data = df
            else:
                data = pd.concat([data, df], ignore_index=True)
        data = data[data["count"] > 0].loc[data["x"] > -1]

        df = None
        for g, d in data.groupby(["feature"]):
            logger.info("Subsampling %s", g)
            x = d["x"].sample(n=sample_size, weights=d["count"], replace=True)
            d = pd.DataFrame({"feature": g[0], "x": x})
            if df is None:
                df = d
            else:
                df = pd.concat([df, d], ignore_index=True)
        output_file = data_path.source_data / output_file_regex.format(
            sampleset=ss
        )
        logger.info("Writing to %s", output_file)
        df[["feature", "x"]].to_csv(
            output_file, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
        )


def _read_nucleotide_diversity(samplesets, features, window_size, regex):
    """Read nucleotide diversity data"""
    if isinstance(window_size, int):
        window_size = [window_size]
    if isinstance(samplesets, str):
        samplesets = [samplesets]

    for ss in samplesets:
        for ws in window_size:
            infile = regex.format(sampleset=ss, window_size=ws)
            logger.info("Reading %s", infile)
            try:
                df = pd.read_csv(infile, sep=",", header=0)
                df = df[df["feature"].isin(features)]
                df["ds"] = ss
                df["window_size"] = ws
                logger.info("Read %s entries", df.shape[0])
            except FileNotFoundError:
                logger.error("File not found: %s", infile)
                continue
            if "data" not in locals():
                data = df
            else:
                data = pd.concat([data, df], ignore_index=True)
    data["feature"] = pd.Categorical(
        data["feature"], categories=features, ordered=True
    )
    data = data.sort_values(
        ["ds", "chrom", "feature", "begin"], ascending=[True, True, True, True]
    )
    return data


@cli.command()
@data_path_opt
@samplesets_opt(default=SAMPLESETS[0:4])
@features_opt(FEATURES)
@window_size_opt("1000000")
@regex_opt(
    "results/diversity/MQ10/{sampleset}/sites.pi.w{window_size}.summary.csv.gz"
)
@paccessible_opt(default=0.0)
@total_sites_opt(default=0)
@verbose_option()
def compile_nucleotide_diversity(
    data_path,
    samplesets,
    features,
    window_size,
    regex,
    paccessible,
    total_sites,
):
    """Compile nucleotide diversity data"""
    logger.info("Compiling nucleotide diversity data")

    data = _read_nucleotide_diversity(samplesets, features, window_size, regex)

    feature_size = get_feature_size(
        data[data["ds"] == samplesets[0]], features
    )
    feature_size = feature_size.sort_values("feature")
    output_file = data_path.source_data / "tbl-feature-size.csv"
    logger.info("Saving feature sizes to %s", output_file)
    feature_size[["feature", "n_sites", "size"]].to_csv(
        output_file, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
    )
    # Make diversity summary table
    results = []
    for g, d in data.groupby(["ds", "feature"]):
        d = d[
            (d["n_accessible"] / d["n_sites"] >= paccessible)
            & (data["n_sites"] >= total_sites)
        ]
        n_total_sites = np.sum(d["n_sites"])
        score = np.sum(d["score"])
        n_accessible = np.sum(d["n_accessible"])

        theta_simple = (
            np.sum(d["n_segregating_sites_accessible"])
            / n_accessible
            / HARMONIC_NUMBER[THETA_CALL_DEPTH[g[0]]]
        )

        results.append(
            pd.Series(
                {
                    "ds": g[0],
                    "feature": g[1],
                    "total_sites": n_total_sites,
                    "total_score": np.sum(d["total_score"]),
                    "score": score,
                    "n_accessible": n_accessible,
                    "p_accessible": n_accessible / n_total_sites,
                    "n_S": np.sum(d["n_segregating_sites_accessible"]),
                    "site_score": score / n_accessible,
                    "theta": np.mean(d["theta"]),
                    "theta_simple": theta_simple,
                }
            )
        )

    nucdiv = pd.DataFrame(results)
    nucdiv["feature"] = pd.Categorical(
        nucdiv["feature"], categories=features, ordered=True
    )
    nucdiv = nucdiv.sort_values(["ds", "feature"], ascending=[True, True])
    output_file = data_path.source_data / "tbl-diversity-summary.csv"
    logger.info("Saving nucleotide diversity summary to %s", output_file)
    nucdiv.to_csv(
        output_file, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
    )
    output_file = data_path.source_data / "tbl-diversity-summary_README.md"
    logger.info("Saving README to %s", output_file)
    with open(output_file, "w") as f:
        info_d = {
            "ds": "Data set",
            "feature": "Feature",
            "total_sites": "Total number of sites for a given feature",
            "total_score": (
                "Aggregate nucleotide diversity score (pi) over all sites"
            ),
            "score": (
                "Aggregate nucleotide diversity score (pi) "
                "over accessible sites"
            ),
            "n_accessible": (
                "Number of accessible sites as determined by coverage masks"
            ),
            "p_accessible": "Proportion of accessible sites",
            "n_S": "Number of accessible segregating sites",
            "site_score": (
                "Mean nucleotide diversity (pi) over all accessible sites"
            ),
            "theta": (
                "Mean Watterson's theta over all accessible sites in windows"
            ),
            "theta_simple": (
                "Watterson's theta calculated as the total number "
                "of accessible sites corrected by average call depth"
            ),
        }
        info = pd.DataFrame.from_dict(info_d, orient="index").reset_index()
        info.columns = ["Field", "Description"]
        f.write(
            "# tbl-diversity-summary.csv\n\n"
            "## Summary of nucleotide diversity statistics\n\n"
        )
        f.write(info.to_markdown(index=False))

    output_file = (
        data_path.source_data / "tbl-pabies-diversity-summary-twocol.csv"
    )
    logger.info(
        "Saving nucleotide diversity two column summary, P.abies, to %s",
        output_file,
    )
    pabies = nucdiv[nucdiv["ds"] == "P.abies"]
    pabies = pabies[["feature", "site_score"]]
    pabies.to_csv(
        output_file, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
    )

    output_file = (
        data_path.source_data / "tbl-highcov-diversity-summary-twocol.csv"
    )
    logger.info(
        "Saving nucleotide diversity two column summary, highcov, to %s",
        output_file,
    )
    highcov = nucdiv[nucdiv["ds"] == "highcov"]
    highcov = pabies[["feature", "site_score"]]
    highcov.to_csv(
        output_file, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
    )


@cli.command()
@output_file_regex_opt(
    "nucleotide-diversity-{sampleset}-w{window_size}.tsv.gz"
)
@data_path_opt
@samplesets_opt(default=SAMPLESETS[0])
@features_opt(FEATURES)
@window_size_opt("100000")
@regex_opt(
    "results/diversity/MQ10/{sampleset}/sites.pi.w{window_size}.summary.csv.gz"
)
@verbose_option()
def compile_genome_scale_nucleotide_diversity(
    data_path, samplesets, features, window_size, regex, output_file_regex
):
    """Compile genome scale nucleotide diversity"""
    logger.info(
        (
            "Compiling genome scale nucleotide diversity data, "
            "all unfiltered (raw) data"
        )
    )

    data = _read_nucleotide_diversity(samplesets, features, window_size, regex)

    outdir = data_path.figshare / "genome-scale-variation"
    if not outdir.exists():
        outdir.mkdir()
    outfile_name = outdir / output_file_regex.format(
        sampleset=samplesets, window_size=window_size
    )

    logger.info("Writing to %s", outfile_name)
    data.to_csv(
        outfile_name, sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC
    )


@cli.command()
@output_file_regex_opt(
    "fig-nucleotide-diversity-{sampleset}-w{window_size}-{paccessible}-{total_sites}.tsv.gz"
)
@data_path_opt
@samplesets_opt(default=SAMPLESETS[0])
@features_opt(FEATURES)
@window_size_opt("100000")
@regex_opt(
    "results/diversity/MQ10/{sampleset}/sites.pi.w{window_size}.summary.csv.gz"
)
@click.option(
    "--paccessible",
    type=float,
    default=0.5,
    help="Minimum proportion of accessible sites",
)
@click.option(
    "--total-sites",
    type=int,
    default=1000,
    help="Minimum number of total sites in a window",
)
@verbose_option()
def compile_filtered_nucleotide_diversity(
    data_path,
    samplesets,
    features,
    window_size,
    regex,
    output_file_regex,
    paccessible,
    total_sites,
):
    """Compile filtered nucleotide diversity for figure"""
    logger.info(
        (
            "Compiling genome scale nucleotide diversity data, "
            "all unfiltered (raw) data"
        )
    )

    data = _read_nucleotide_diversity(samplesets, features, window_size, regex)
    data["window_size"] = data["end"] - data["begin"]
    # Filter data
    data = data[
        (data["n_accessible"] / data["n_sites"] >= paccessible)
        & (data["n_sites"] >= total_sites)
    ]

    outfile_name = data_path.source_data / output_file_regex.format(
        sampleset=samplesets,
        window_size=window_size,
        paccessible=paccessible,
        total_sites=total_sites,
    )

    logger.info("Writing to %s", outfile_name)
    data = data[
        ["site_score", "feature", "n_accessible", "n_sites", "window_size"]
    ]
    data.to_csv(
        outfile_name, sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC
    )


@cli.command()
@output_file_regex_opt("fig-diversity-genic-proximity-{sampleset}.csv.gz")
@data_path_opt
@samplesets_opt(default=SAMPLESETS[1])
@features_opt(["CDS", "gene"])
@window_size_opt("1000")
@regex_opt(
    "results/diversity/MQ10/{sampleset}/sites.pi.w{window_size}.e200000.{feature}.summary.csv.gz"
)
@click.option(
    "--paccessible",
    type=float,
    default=0.5,
    help="Minimum proportion of accessible sites",
)
@verbose_option()
def compile_genic_proximity(  # pylint: disable=too-many-locals
    data_path,
    samplesets,
    features,
    window_size,
    regex,
    paccessible,
    output_file_regex,
):
    """Compile nucleotide diversity in genic proximity for figure"""
    logger.info("Compiling nucleotide diversity genic proximity data")

    def _wpos(x):
        fac = 1 if x[0] == "D" else -1
        return fac * int(x[1:]) + fac * 0.5

    def _bootstrap(data, size=200):
        x = bootstrap((data,), np.mean, n_resamples=size, method="percentile")
        result = pd.Series(
            [
                np.mean(data),
                np.mean(x.bootstrap_distribution),
                x.confidence_interval.low,
                x.confidence_interval.high,
            ],
            index=["mean", "boot.mean", "ci.low", "ci.high"],
        )
        return result

    infile_regex = [
        regex.format(
            sampleset="{sampleset}", window_size="{window_size}", feature=ft
        )
        for ft in features
    ]
    for i, ft in enumerate(features):
        df = _read_nucleotide_diversity(
            samplesets, ["genome"], window_size, infile_regex[i]
        )
        df["feature"] = ft
        df["window"] = df["region"].str.replace(
            r".+_([UD]\d+)$", r"\1", regex=True
        )
        df["geneid"] = (
            df["region"]
            .str.replace(r"_[UD]\d+$", r"", regex=True)
            .astype("category")
            .cat.codes.map(lambda x: f"G{x}")
        )
        df["wpos"] = df["window"].map(_wpos)
        df = df[(df["n_accessible"] / df["n_sites"] >= paccessible)]

        logger.info("Calculating bootstrap statistics for pi, feature %s", ft)
        bsdata1 = (
            df.groupby(["wpos", "feature"])["site_score"]
            .apply(_bootstrap)
            .rename_axis(index=["wpos", "feature", "stat"])
            .reset_index()
            .pivot_table(
                index=["wpos", "feature"], columns="stat", values="site_score"
            )
        )
        bsdata1["statistic"] = "pi"
        logger.info(
            "Calculating bootstrap statistics for theta, feature %s", ft
        )
        bsdata2 = (
            df.groupby(["wpos", "feature"])["theta"]
            .apply(_bootstrap)
            .rename_axis(index=["wpos", "feature", "stat"])
            .reset_index()
            .pivot_table(
                index=["wpos", "feature"], columns="stat", values="theta"
            )
        )
        bsdata2["statistic"] = "theta"
        bsdata = pd.concat([bsdata1, bsdata2]).reset_index()

        if "data" not in locals():
            data = bsdata
        else:
            data = pd.concat([data, bsdata])

    outfile_name = data_path.source_data / output_file_regex.format(
        sampleset=samplesets,
        window_size=window_size,
    )
    logger.info("Writing to %s", outfile_name)
    data.to_csv(
        outfile_name, sep=",", index=False, quoting=csv.QUOTE_NONNUMERIC
    )
