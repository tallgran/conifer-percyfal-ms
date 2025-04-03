"""Convert FASTAINDEX to BEDFILE format."""

import logging

import click

logger = logging.getLogger(__name__)


@click.command(help=__doc__)
@click.argument("fastaindex", type=click.File("r"))
@click.argument("bedfile", type=click.File("w"))
def cli(fastaindex, bedfile):
    while True:
        chunk = fastaindex.readline()
        if not chunk:
            break
        data = chunk.split("\t")
        bedfile.write(f"{data[0]}\t0\t{data[1]}\n")
