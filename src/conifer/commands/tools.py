"""Conifer tools and scripts.

Tools that also can be called as standalone scripts by adding prefix
'conifer-' to name in list below.

"""

import logging
import os

import click

from conifer.cli import CONTEXT_SETTINGS, ConiferCli

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]


logger = logging.getLogger(__name__)


class ConiferToolsCli(ConiferCli):
    """Conifer tools command line interface"""

    cmd_folder = os.path.abspath(
        os.path.join(os.path.dirname(__file__), os.pardir, "tools")
    )
    module = "conifer.tools"


@click.command(
    cls=ConiferToolsCli,
    context_settings=CONTEXT_SETTINGS,
    help=__doc__,
    name="tools",
)
def cli():
    """Conifer tools and scripts"""
    logger.info(__shortname__)
