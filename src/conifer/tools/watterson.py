"""Watterson

Calculate Watterson's theta, either directly from a VCF file or from a
summary file.

Corrections:

- average number of samples -> corrected harmonic number

"""

import logging
import sys

import click
import numpy as np
import pandas as pd

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)

HARMONIC_NUMBER_MAX = 2500

# Harmonic number. For large n, the naive function is slow and an
# approximation should be used (cf
# https://stackoverflow.com/questions/404346/python-program-to-calculate-harmonic-series).
# The largest n we see is <2100 so we can calculate a lookup
# dictionary to speed up calculations.
HARMONIC_NUMBER = np.zeros(HARMONIC_NUMBER_MAX + 1)
HARMONIC_NUMBER_SQUARED = np.zeros(HARMONIC_NUMBER_MAX + 1)
HARMONIC_NUMBER[0] = np.nan
HARMONIC_NUMBER_SQUARED[0] = np.nan
for ii in np.arange(1, HARMONIC_NUMBER_MAX):
    HARMONIC_NUMBER[ii] = np.sum(1 / np.arange(1, ii))
    HARMONIC_NUMBER_SQUARED[ii] = np.sum(1 / (np.arange(1, ii) ** 2))

EULER_MASCHERONI = 0.57721566490153286060651209008240243104215933593992


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """Utilities for calculating Watterson's theta"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.debug("Running %s version %s", __shortname__, __version__)


@cli.command()
@click.argument("snpeff_table", type=click.Path(exists=True))
@click.argument("kaks_table", type=click.Path(exists=True))
@click.option(
    "--n-samples",
    "-n",
    help="number of samples",
    default=(1000,),
    type=int,
    multiple=True,
)
@click.option(
    "--proportion-accessible",
    "-p",
    help="proportion accessible sites",
    default=1.0,
    type=click.FloatRange(0.0, 1.0),
)
@click.option("--dataset", "-d", help="dataset label for output", default="ds")
@click.option(
    "show_header", "--header", "-H", help="header", default=False, is_flag=True
)
@click.version_option(version=__version__)
@click.pass_context
def snpeff(
    ctx,
    snpeff_table,
    kaks_table,
    proportion_accessible,
    n_samples,
    dataset,
    show_header,
):
    """Calculate Watterson's theta from snpEff summary file."""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.debug("Running %s version %s", __shortname__, __version__)

    # Read snpEff summary file
    with open(snpeff_table) as f:
        header = f.readline().strip().split("\t")
    if header[0] != "#GeneName":
        data = pd.read_table(snpeff_table, sep="\t", skiprows=1)
    else:
        data = pd.read_table(snpeff_table, sep="\t")

    data.rename(columns={"#GeneName": "GeneName"}, inplace=True)
    data = data[
        [
            "TranscriptId",
            "variants_effect_missense_variant",
            "variants_effect_synonymous_variant",
        ]
    ]

    # Read KaKs summary file
    kaks = pd.read_table(kaks_table, sep="\t")
    kaks = kaks[["Sequence", "S-Sites", "N-Sites"]]
    kaks.rename(columns={"Sequence": "TranscriptId"}, inplace=True)

    # Merge the two tables
    merged = pd.merge(data, kaks, on="TranscriptId")
    merged.to_csv("merged.csvtk.csv", index=False)
    res = []

    # Calculate Watterson's theta
    for n in n_samples:
        h_n = (
            HARMONIC_NUMBER[n]
            if n < HARMONIC_NUMBER_MAX
            else np.log(n) + EULER_MASCHERONI
        )

        n_s = np.sum(merged["variants_effect_synonymous_variant"])
        n_n = np.sum(merged["variants_effect_missense_variant"])
        n_ns = n_s + n_n

        s_tot = np.sum(merged["S-Sites"])
        n_tot = np.sum(merged["N-Sites"])
        ns_tot = s_tot + s_tot

        theta_s_raw = n_s / s_tot / h_n
        theta_n_raw = n_n / n_tot / h_n
        theta_ns_raw = n_ns / ns_tot / h_n

        theta_s = theta_s_raw / proportion_accessible
        theta_n = theta_n_raw / proportion_accessible
        theta_ns = theta_ns_raw / proportion_accessible

        res.append(
            (
                dataset,
                n,
                h_n,
                proportion_accessible,
                n_s,
                n_n,
                n_ns,
                s_tot,
                n_tot,
                ns_tot,
                theta_s_raw,
                theta_n_raw,
                theta_ns_raw,
                theta_s,
                theta_n,
                theta_ns,
            )
        )

    x = pd.DataFrame(
        res,
        columns=[
            "dataset",
            "n",
            "H(n)",
            "pacc",
            "nS",
            "nN",
            "nNS",
            "Stot",
            "Ntot",
            "NStot",
            "thetaS_raw",
            "thetaN_raw",
            "thetaNS_raw",
            "thetaS",
            "thetaN",
            "thetaNS",
        ],
    )
    x["n"] = x["n"].astype(int)
    if show_header:
        x.to_csv(sys.stdout, index=False, sep="\t")
    else:
        x.to_csv(sys.stdout, index=False, header=False, sep="\t")


@cli.command()
@click.argument("vcf_summary", type=click.Path(exists=True))
@click.option(
    "--n-samples",
    "-n",
    help="number of samples",
    default=(1000,),
    type=int,
    multiple=True,
)
@click.version_option(version=__version__)
@click.pass_context
def vcf_summary(ctx, vcf_summary, n_samples):
    """Calculate Watterson's theta from a VCF summary file."""
    data = pd.read_table(vcf_summary, sep=",")
    gcols = ["feature", "statistic"]
    sitecols = [
        "score",
        "total_score",
        "n_sites",
        "n_accessible",
        "n_segregating_sites",
        "n_segregating_sites_accessible",
    ]
    divcols = ["site_score", "theta"]

    print("Site statistics, averaged over regions, no filtering")
    df = data[gcols + divcols].groupby(by=gcols).mean()
    print(df)

    print("Site statistics, summed over regions, no filtering")
    df = data[gcols + sitecols].groupby(by=gcols).sum()
    for n in n_samples:
        h_n = (
            HARMONIC_NUMBER[n]
            if n < HARMONIC_NUMBER_MAX
            else np.log(n) + EULER_MASCHERONI
        )
        df[f"theta_{n}"] = (
            df["n_segregating_sites_accessible"] / df["n_accessible"] / h_n
        )
    print(df)


@cli.command()
@click.argument("vcf_summary", type=click.Path(exists=True))
@click.option(
    "--n-samples",
    "-n",
    help="number of samples",
    default=(1000,),
    type=int,
    multiple=True,
)
@click.version_option(version=__version__)
@click.pass_context
def vcf_summary_qc(ctx, vcf_summary, n_samples):
    """Perform QC on vcf summary file"""
    data = pd.read_table(vcf_summary, sep=",")
    sitecols = [
        "score",
        "total_score",
        "n_sites",
        "n_accessible",
        "n_segregating_sites",
        "n_segregating_sites_accessible",
    ]

    for n in n_samples:
        h_n = (
            HARMONIC_NUMBER[n]
            if n < HARMONIC_NUMBER_MAX
            else np.log(n) + EULER_MASCHERONI
        )
        data[f"theta_qc_{n}"] = (
            data["n_segregating_sites_accessible"] / data["n_accessible"] / h_n
        )
        x = abs((data[f"theta_qc_{n}"] - data["theta"]) / data["theta"])
        i = sum(x > 0.05)
        print(f"Number of sites with >5% difference in theta_{n}: {i}")
        print(data[x > 0.05][sitecols])
