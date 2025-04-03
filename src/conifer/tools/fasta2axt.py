"""fasta2axt

Convert fasta file to axt format.
"""

import logging

import click
import pyfastx

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)


@click.command(help=__doc__, name=__shortname__)
@click.argument("fasta", type=click.Path(exists=True))
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx, fasta):
    """Convert a FASTA file to AXT."""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__
    logger.debug("Running %s version %s", __shortname__, __version__)

    fa = pyfastx.Fasta(fasta)
    for seq in fa:
        print(seq.name)
        print(seq.seq)
        print(seq.seq)
        print()
