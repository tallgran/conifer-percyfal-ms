"""Generate data for manuscript submission.

Create the source data and figshare data for the manuscript
submission.

"""

import logging

import click
from nbis import wrappers
from nbis.cli import pass_environment
from nbis.snakemake import (
    jobs_option,
    no_profile_option,
    profile_option,
    snakemake_argument_list,
)

from conifer.snakemake import config

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]


logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
def cli():
    """Main CLI entry point."""
    logger.debug("Running %s", __shortname__)


@cli.command(
    context_settings={"ignore_unknown_options": True},
    help="Generate manuscript data",
)
@profile_option()
@no_profile_option()
@jobs_option()
@snakemake_argument_list()
@pass_environment
def run(env, profile, no_profile, jobs, snakemake_args):
    """Generate manuscript data"""
    smk_options = list(snakemake_args) + jobs + profile
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "manuscript-run.smk"
    wrappers.snakemake(options=smk_options, snakefile=snakefile, targets="")


@cli.command(context_settings={"ignore_unknown_options": True}, help="qc")
@profile_option(default="local")
@no_profile_option()
@jobs_option()
@snakemake_argument_list()
@pass_environment
def qc(env, profile, no_profile, jobs, snakemake_args):
    """Command docstring"""
    smk_options = list(snakemake_args) + jobs + profile
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "manuscript-qc.smk"
    wrappers.snakemake(
        options=" ".join(smk_options), snakefile=snakefile, targets=""
    )


@cli.command(
    context_settings={"ignore_unknown_options": True}, help="eva-submission"
)
@profile_option(default="local")
@no_profile_option()
@jobs_option()
@snakemake_argument_list()
@pass_environment
def eva_submission(env, profile, no_profile, jobs, snakemake_args):
    """"""
    smk_options = list(snakemake_args) + jobs + profile
    snakefile = (
        config.SNAKEMAKE_ROOT / "commands" / "manuscript-eva-submission.smk"
    )
    wrappers.snakemake(
        options=" ".join(smk_options), snakefile=snakefile, targets=""
    )
