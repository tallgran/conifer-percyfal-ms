"""Console script for conifer

Collection of scripts and workflows to run conifer analyses.

## Installation

See README.md for installation instructions.

## Quick usage

The conifer toolkit is run from the command line using the conifer
command, which will provide a list of available subcommands. Prior to
any analysis, a number of configuration files need to be created. The
administrative commands are:

- conifer datasources - setup external datasources via links/rsync.
  See datasources config file for details.

## Tools

The tools subcommand serves as an entry point to standalone tools and
scripts. They can also be accessed as separate commands by prefixing
the subcommand with 'conifer-', e.g., conifer-summarize-diversity.

## Smk

The smk subcommand provides a list of snakemake workflows which are to
be run in the following order:

0a. conifer smk init-samplesets - initialize sample sets

0b. conifer smk annotation - generate region definition files
    (feature) for diversity analyses

1. conifer smk coverage - run coverage analyses for all samples

2. conifer smk agregate-coverage - aggregate sample coverages for
   different sample sets

3. conifer smk diversity - create sequence masks based on output from
   step 2. and calculate diversity statistics for sample sets

### Additional analyses

1. conifer smk watterson - calculate Watterson's theta for synonymous
and non-synonymous sites using SnpEff annotations and KaKs

2. conifer smk presabs - calculate presence/absence of genes in sample
sets to prepare for notebook/presabs.ipynb
"""

import logging
import os
import pathlib

import click
from nbis import decorators
from nbis.config import load_config
from nbis.env import Environment

from . import __version__

__author__ = "Per Unneberg"

logger = logging.getLogger(__name__)


PKG_DIR = pathlib.Path(__file__).absolute().parent
CONTEXT_SETTINGS = {"auto_envvar_prefix": "CONIFER", "show_default": True}

pass_environment = click.make_pass_decorator(Environment, ensure=True)


class ConiferCli(click.MultiCommand):
    """Conifer command line interface"""

    module = "conifer.commands"
    cmd_folder = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "commands")
    )

    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir(self.cmd_folder):
            if (
                filename.endswith(".py")
                and not filename.startswith("__")
                and not filename.startswith(".")
            ):
                rv.append(filename[:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, cmd_name):
        mod = __import__(f"{self.module}.{cmd_name}", None, None, ["cli"])
        return mod.cli


@click.command(
    cls=ConiferCli,
    context_settings=CONTEXT_SETTINGS,
    help=__doc__,
    name="conifer",
)
@click.version_option(version=__version__)
@click.option(
    "--config-file", help="configuration file", type=click.Path(exists=True)
)
@decorators.debug_option()
@pass_environment
def cli(env, config_file):
    """Conifer command line interface"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
    )
    if env.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if config_file is None:
        env.home = pathlib.Path(os.curdir).absolute()
        config_file = env.home / "conifer.yaml"
    else:
        config_file = pathlib.Path(config_file).absolute()
        env.home = config_file.parent
    if config_file.exists():
        config = load_config(file=config_file)
    else:
        config = load_config(data={"project_name": "conifer"})
    env.config = config
