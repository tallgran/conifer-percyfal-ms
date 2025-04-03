"""Run snakemake workflows.

Main entry point to run snakemake workflows.
"""

import glob
import logging
import os
import shutil
import subprocess as sp
from pathlib import Path

import click
from nbis import wrappers
from nbis.snakemake import (
    directory_option,
    envmodules_configfile_option,
    jobs_option,
    no_profile_option,
    profile_option,
    report_option,
    snakemake_argument_list,
)

from conifer.cli import pass_environment
from conifer.snakemake import config

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]


logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
def cli():
    """Run snakemake workflows"""
    logger.debug("Running %s", __shortname__)


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option()
@no_profile_option()
@jobs_option()
@report_option(report_file=Path("reports/coverage/report.html"))
@directory_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def coverage(
    env,
    profile,
    no_profile,
    jobs,
    report,
    directory,
    envmodules_configfile,
    snakemake_args,
):  # pylint: disable=R0913
    """Run coverage analysis.

    The goal of this analysis is to generate coverage files and benchmark
    programs and file formats. The command runs both snakemake analyses
    and a regular script.
    """
    webexport = os.path.join(env.config.webexport.url, "qc/coverage")
    options = (
        list(snakemake_args)
        + jobs
        + profile
        + report
        + directory
        + envmodules_configfile
    )
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "coverage.smk"

    wrappers.snakemake(options=options, snakefile=snakefile, targets="")
    if len(report) > 0:
        shutil.copy2(report[1], webexport)
        report_dir = Path("reports/coverage")
        mqc = glob.glob(str(report_dir / "multiqc*"))
        sp.run(["rsync", "-av"] + mqc + [webexport], check=True)


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option()
@no_profile_option()
@jobs_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@report_option(report_file="reports/annotation/report.html")
@snakemake_argument_list()
@pass_environment
def annotation(
    env,
    profile,
    no_profile,
    jobs,
    envmodules_configfile,
    report,
    snakemake_args,
):
    """Convert annotation files to regions in data/resources/regions.

    Convert annotation files to regions. Outputs are found in
    data/resources/regions.

    Annotation coordinates vs mapping coordinates

    Original annotation coordinates are defined wrt
    Picab02_chromosomes.fasta.gz. However, the coordinates of mapped
    sequences and consequently variants are defined wrt pabies-2.0.fa,
    which is the locked reference sequence but with chopped up
    genomes. Therefore, annotation coordinates have been reassigned in
    files named *.liftover.gff3.gz and similar.

    """
    options = (
        list(snakemake_args) + jobs + profile + envmodules_configfile + report
    )
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "smk-annotation.smk"
    wrappers.snakemake(options=options, snakefile=snakefile, targets="")


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option()
@no_profile_option()
@jobs_option()
@report_option(report_file=Path("reports/diversity/report.html"))
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def diversity(
    env,
    profile,
    no_profile,
    jobs,
    envmodules_configfile,
    snakemake_args,
    report,
):
    """Calculate genetic diversity with vcftools and friends.

    Calculate genetic diversity (pi, theta) and run selection scans
    (Tajima's D) based on sites defined by mask files from mask
    module.

    Window summaries are compiled with conifer-summarize-diversity.
    """
    webexport = os.path.join(env.config.webexport.url, "qc/diversity")
    options = (
        list(snakemake_args) + jobs + profile + report + envmodules_configfile
    )
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "smk-diversity.smk"
    wrappers.snakemake(options=options, snakefile=snakefile, targets="")
    if len(report) > 0:
        shutil.copy2(report[1], webexport)


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option()
@no_profile_option()
@jobs_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def qc(env, profile, no_profile, jobs, envmodules_configfile, snakemake_args):
    """Collect qc metrics"""
    options = list(snakemake_args) + jobs + profile + envmodules_configfile
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "smk-qc.smk"
    wrappers.snakemake(options=options, snakefile=snakefile, targets="")


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option()
@no_profile_option()
@jobs_option()
@report_option(report_file="reports/aggregate_coverage/report.html")
@directory_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def aggregate_coverage(
    env,
    profile,
    no_profile,
    jobs,
    directory,
    envmodules_configfile,
    snakemake_args,
    report,
):
    """Generate aggregate coverage files and compute summary statistics.

    Aggregate individual coverage files, grouped by sample sets. Two
    aggregate values are calculated for each site:

    count : count the number of samples that have a minimum
    coverage given by the tag ge{MINIMUM}

    sum : aggregate the total coverage over all samples

    ## Coverage masks

    Coverage masks are generated as the intersection of two sources.

    1. Overall coverage

    Output:
    data/mosdepth_coverage/MQ{MQ}/{sampleset}/
        sum_coverage.filter_l{lower}-u{upper}.p{100}.bed.gz.
    Sites with total coverage outside range lower-upper are filtered out.

    2. Number of samples with minimum coverage

    Output:
    data/mosdepth_coverage/MQ{MQ}/{sampleset}/
        count_l{minalleles}_coverage.filter_l{lower}-u{upper}.p{100}.bed.gz.
    Each site is a count of the number of individuals with a coverage >=
    minalleles. Sites with at least `lower` number of individuals are kept.

    The output is put in data/mask/MQ{MQ}/{sampleset} and consists of

    1. an interval file with kept intervals

    sum_coverage.filter_l4248-u11258--count_l3_coverage.filter_l528-u1070.bed.gz

    2. a fasta mask file

    sum_coverage.filter_l4248-u11258--count_l3_coverage.filter_l528-u1070.single.fa

    which is used as input to vcftools.


    ## Region masks

    Region masks consist of intervals from annotation files located in
    data/resources:

    - Picab02_codingAll.gff3.gz
    - Picab02_repeats.gff3.gz
    - Picab02_TEs.gff3.gz

    Parsed regions are in data/resources/regions.

    Parsed regions are combined with coverage masks to output files with
    prefix
    {region}.sum_coverage.filter_l{lower1}-u{upper1}--
        count_l{min_alleles}_coverage.filter_l{lower2}-u{upper2}
    """
    webexport = os.path.join(env.config.webexport.url, "qc/aggregate_coverage")
    options = (
        list(snakemake_args)
        + jobs
        + profile
        + report
        + directory
        + envmodules_configfile
    )
    snakefile = (
        config.SNAKEMAKE_ROOT / "commands" / "smk-aggregate-coverage.smk"
    )
    wrappers.snakemake(
        options=" ".join(options), snakefile=snakefile, targets=""
    )
    if len(report) > 0:
        shutil.copy2(report[1], webexport)


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option()
@no_profile_option()
@jobs_option()
@directory_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def init_samplesets(
    env,
    profile,
    no_profile,
    jobs,
    directory,
    envmodules_configfile,
    snakemake_args,
):
    """Initialize sample sets.

    Initialize sample sets by creating a samples.tsv file in
    resources/samplesets/{sampleset}.

    """
    options = (
        list(snakemake_args)
        + jobs
        + profile
        + directory
        + envmodules_configfile
    )
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "smk-init-samplesets.smk"
    wrappers.snakemake(options=options, snakefile=snakefile, targets="")


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option(default="local")
@no_profile_option()
@jobs_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def presabs(
    env, profile, no_profile, jobs, snakemake_args, envmodules_configfile
):
    """Run presence/absence workflow"""
    smk_options = list(snakemake_args) + jobs + profile + envmodules_configfile
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "smk-presabs.smk"
    wrappers.snakemake(
        options=" ".join(smk_options), snakefile=snakefile, targets=""
    )


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option(default="local")
@no_profile_option()
@jobs_option()
@envmodules_configfile_option(envmodules_configfile="config/envmodules.yaml")
@snakemake_argument_list()
@pass_environment
def watterson(
    env, profile, no_profile, jobs, snakemake_args, envmodules_configfile
):
    """Run Watterson workflow"""
    smk_options = list(snakemake_args) + jobs + profile + envmodules_configfile
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "smk-watterson.smk"
    wrappers.snakemake(
        options=" ".join(smk_options), snakefile=snakefile, targets=""
    )
