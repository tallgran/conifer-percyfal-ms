"""Run datasources setup.

Run the datasources snakemake workflow to setup datasources for
project. The workflow expects a 'config/datasources.yaml'
configuration file.
"""

import logging

import click
from nbis import wrappers
from nbis.cli import pass_environment
from nbis.snakemake import jobs_option, no_profile_option, profile_option

from conifer.snakemake import config

logger = logging.getLogger(__name__)


@click.group(
    name="datasources",
    help=__doc__,
)
@click.pass_context
def cli(ctx):
    """Run datasources setup."""
    logger.debug(f"Running {ctx.find_root().info_name} {ctx.info_name}")


@cli.command(context_settings={"ignore_unknown_options": True})
@profile_option(default="local")
@no_profile_option()
@jobs_option()
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
@pass_environment
def run(env, profile, no_profile, jobs, snakemake_args):
    """Setup datasources using snakemake v8 storage plugins."""
    smk_options = list(snakemake_args) + jobs + profile
    snakefile = config.SNAKEMAKE_ROOT / "commands" / "datasources-run.smk"
    wrappers.snakemake(
        options=" ".join(smk_options), snakefile=snakefile, targets=""
    )
