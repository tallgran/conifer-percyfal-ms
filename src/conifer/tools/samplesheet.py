"""Conifer samplesheet.

Make samplesheet from source defined in environment variable
SAMPLESHEET_URL. Output is written to file defined by variable SAMPLE_INFO_OUT.


Data was downloaded as csv file and synced to uppmax:

\b
rsync -av /home/peru/Downloads/SAMPLE_INFO\\ -\\ Spruce-samples.csv
    /proj/uppstore2017145/V3/users/perun/conifer/resources/

The script will add an additional column grouping samples into
breeding population, range wide samples, and outgroup samples.

20221025: four missing samples mapped but as of yet not variant
called. Outgroup VCF-ID names P_engelmannii, P_glauca, and
P_sitchensis do not match names in vcf files, Engelmannii_SRX1530215,
Glauca_SRX160982, and Sitchensis_SRX2354260.

\b
Obsolete
========


FIXED 20221025: 20220916: Samples on rows 1086-1095 are misaligned and
have to be "caught" with regex
KAM|KRY|SOC|GUL|SMU|ORM|TAR|PEL|BOR|ACR_16H311 and have to be mapped
to sample ids P21002-133 - P21002-142. They lack metadata and need
fixing.

FIXED 20221025: Samples P21002_116 and P22451-251 lack species
information. Assuming Picea abies for now.

PARTIALLY FIXED 20221025: Outgroup samples Engelmannii_SRX1530215,
Glauca_SRX160982, and Sitchensis_SRX2354260 are also missing.
20221025: outgroup samples present but named P_engelmannii, P_glauca,
and P_sitchensis

FIXED 20221025: Samples P15766-177, P15766-222, P15766-281, and
P21002-128 are unmapped and are removed at this stage to facilitate
downstream processing.

"""

import logging
import pathlib
import re

import click
import pandas as pd

from conifer.cli import pass_environment

logger = logging.getLogger(__name__)


OUTGROUP_SAMPLES = [
    "Engelmannii_SRX1530215",
    "Glauca_SRX160982",
    "Sitchensis_SRX2354260",
    "P.obovata_P24355",
    "P267-101",
]
REFERENCE_SAMPLES = [
    "diploid",
    "haploid_ERX242653",
]
REFERENCE_SAMPLES_FROM = ["Diploid", "Haploid_ERX242654"]
OUTGROUP_SAMPLES_FROM = [
    "P_engelmannii",
    "P_glauca",
    "P_sitchensis",
    "P24355-101",
]
OUTGROUP_SAMPLES_TO = [
    "Engelmannii_SRX1530215",
    "Glauca_SRX160982",
    "Sitchensis_SRX2354260",
    "P.obovata_P24355",
]
outgroups = "|".join(OUTGROUP_SAMPLES_FROM)
refs = "|".join(REFERENCE_SAMPLES_FROM)
SAMPLE_RE = re.compile(
    rf"^(P[0-9]+[\-_][0-9]+|P.obovata_P24355|{outgroups}|{refs})$"
)
SAMPLE_INFO = pathlib.Path("resources") / "SAMPLE_INFO - Spruce-samples.csv"
SAMPLE_INFO_OUT = pathlib.Path("resources") / "samplesheet.tsv"

if not SAMPLE_INFO.exists():
    logger.error("Could not find '%s'", SAMPLE_INFO)
    logger.error("Please download from source and sync to HPC")
    raise FileNotFoundError("Could not find '%s'" % SAMPLE_INFO)


@click.command(help=__doc__)
@click.option("--outfile", "-o", type=click.File("w"), default=SAMPLE_INFO_OUT)
@click.option("--show", is_flag=True, help="print samplesheet to screen")
@pass_environment
def cli(env, outfile, show):
    click.echo(f"Reading samplesheet from {SAMPLE_INFO}")
    if not show:
        click.echo(f"Writing samplesheet to {outfile}")
    else:
        click.echo("Showing samplesheet on screen")
    info = pd.read_csv(env.home / SAMPLE_INFO)
    columns = [
        "Scilife-ID",
        "VCF-ID",
        "Internal (Lab)-ID",
        "ENA project id",
        "Experiment ID",
        "Biosample",
        "Species",
        "Geographic origin (country)",
        "Population",
        "Lattitude",
        "Longitude",
        "Altitude",
        "Lab",
        "Note",
    ]
    ix = info["VCF-ID"].map(lambda x: SAMPLE_RE.match(x) is not None)
    samples = info.loc[ix][columns]
    samples.insert(loc=1, column="SM", value=samples["VCF-ID"])
    samples.columns = [
        "Scilife-ID",
        "SM",
        "VCF-ID",
        "Lab-ID",
        "ENA-project-ID",
        "Experiment-ID",
        "Biosample",
        "Species",
        "Country",
        "Population",
        "Latitude",
        "Longitude",
        "Altitude",
        "Lab",
        "Note",
    ]
    samples = samples.fillna("x")
    samples = samples.replace("x", None)

    # Rename some samples for consistency
    samples.loc[samples.SM == OUTGROUP_SAMPLES_FROM[0], "SM"] = (
        OUTGROUP_SAMPLES_TO[0]
    )
    samples.loc[samples.SM == OUTGROUP_SAMPLES_FROM[1], "SM"] = (
        OUTGROUP_SAMPLES_TO[1]
    )
    samples.loc[samples.SM == OUTGROUP_SAMPLES_FROM[2], "SM"] = (
        OUTGROUP_SAMPLES_TO[2]
    )
    samples.loc[samples.SM == OUTGROUP_SAMPLES_FROM[3], "SM"] = (
        OUTGROUP_SAMPLES_TO[3]
    )
    samples.loc[samples.SM == REFERENCE_SAMPLES_FROM[0], "SM"] = (
        REFERENCE_SAMPLES[0]
    )
    samples.loc[samples.SM == REFERENCE_SAMPLES_FROM[1], "SM"] = (
        REFERENCE_SAMPLES[1]
    )

    samples["PopulationType"] = "BreedingPopulation"
    samples.loc[
        samples.Population != "Northern breeding population", "PopulationType"
    ] = "RangeWide"
    samples.loc[samples.SM.isin(OUTGROUP_SAMPLES), "PopulationType"] = (
        "Outgroup"
    )
    samples.loc[samples.SM.isin(REFERENCE_SAMPLES), "PopulationType"] = (
        "Reference"
    )

    if show:
        print(samples)
    else:
        samples.to_csv(outfile, sep="\t", index=False, na_rep="NA")
