"""Generate windows around genes for genic proximity analysis.

Generate windows upstream and downstream of genes based on an
ANNOTATION file. A fasta index formatted GENOME file is required for
setting genome sizes. The output file is a BED-formatted file with
windows named after the gene in question and tagged with a direction
("u" for upstream, "d" for downstream) and window index.

"""

import logging
from dataclasses import dataclass

import click
import numpy as np
import pandas as pd
from tqdm import tqdm

from conifer.io.bed import read_fai_genomefile
from conifer.io.files import is_gzip, sniff_infile

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
)


np.seterr(divide="ignore", invalid="ignore")


@dataclass
class Annotation:
    """Annotation dataclass"""

    identifier: str
    chrom: str
    _begin_coord: int
    _end_coord: int
    strand: str

    @property
    def id(self):
        """Return identifier"""
        return self.identifier

    @property
    def begin(self):
        """Return begin of annotation unit"""
        return self._begin_coord

    @property
    def end(self):
        """Return end of annotation unit"""
        return self._end_coord


@dataclass
class Gene(Annotation):
    """Gene dataclass"""

    transcripts: dict
    mode: str = "cds"

    @property
    def begin(self):
        """Return begin of annotation unit"""
        if self.mode == "cds":
            try:
                return min(
                    y.begin for x in self.transcripts.values() for y in x.cds
                )
            except ValueError as e:
                logger.warning(e)
                logger.warning(
                    (
                        "Gene %s (%s:%i-%i) has no cds regions; "
                        "likely overlaps chunk border"
                    ),
                    self.id,
                    self.chrom,
                    self._begin_coord,
                    self._end_coord,
                )
                return self._begin_coord
        return self._begin_coord

    @property
    def end(self):
        """Return end of annotation unit"""
        if self.mode == "cds":
            try:
                return max(
                    y.end for x in self.transcripts.values() for y in x.cds
                )
            except ValueError as e:
                logger.warning(e)
                logger.warning(
                    (
                        "Gene %s (%s:%i-%i) has no cds regions; "
                        "likely overlaps chunk border"
                    ),
                    self.id,
                    self.chrom,
                    self._begin_coord,
                    self._end_coord,
                )
                return self._end_coord
        return self._end_coord


@dataclass
class Transcript(Annotation):
    """Transcript dataclass"""

    parent: str
    cds: list


@dataclass
class CodingSequence(Annotation):
    """Coding sequence (CDS) dataclass"""

    parent: str


def parse_gff_attributes(attributes_str):
    """Parse GFF attributes into a dictionary."""
    attributes_dict = {}
    attributes_list = attributes_str.split(";")
    for attribute in attributes_list:
        if "=" in attribute:
            key, value = attribute.split("=")
            attributes_dict[key] = value
    return attributes_dict


def process_chromosome(  # pylint: disable=too-many-locals,too-many-branches
    df: pd.DataFrame, chrom_length: int, extension: int, gene_flag: bool
) -> pd.DataFrame:
    """Process chromosome to find CDS (or genes) with at least
    EXTENSION bp apart."""
    genes = {}
    transcripts = {}
    dist = []
    gene = None
    dext = 2 * extension
    for i, row in df.iterrows():
        attr = parse_gff_attributes(row.attribute)
        if row.feature == "gene":
            gene = Gene(
                attr["ID"], row.chrom, row.begin, int(row.end), row.strand, {}
            )
            genes[gene.id] = gene
            if gene_flag:
                gene.mode = "gene"
        elif row.feature == "mRNA":
            mrna = Transcript(
                attr["ID"],
                row.chrom,
                row.begin,
                int(row.end),
                row.strand,
                attr["Parent"],
                [],
            )
            if mrna.parent not in genes:
                continue
            transcripts[mrna.id] = mrna
            genes[mrna.parent].transcripts[mrna.id] = mrna
        elif row.feature == "CDS":
            cds = CodingSequence(
                attr["ID"],
                row.chrom,
                row.begin,
                int(row.end),
                row.strand,
                attr["Parent"],
            )
            if cds.parent not in transcripts:
                continue
            gene_id = transcripts[cds.parent].parent
            genes[gene_id].transcripts[cds.parent].cds.append(cds)
    if len(genes) == 0:
        return None
    genes = list(genes.values())
    if len(genes) == 1:
        dist = [genes[0].begin, chrom_length - genes[0].end]
    else:
        dist = [genes[0].begin]
        for i, gene in enumerate(genes[1:]):
            try:
                dist.append(gene.begin - genes[i - 1].end)
            except ValueError as e:
                logger.warning(e)
                logger.warning(
                    (
                        "Gene %s (%s:%i-%i) has no cds regions; "
                        "likely overlaps chunk border"
                    ),
                    gene.id,
                    gene.chrom,
                    gene.begin,
                    gene.end,
                )
        try:
            dist.append(chrom_length - genes[-1].end)
        except ValueError as e:
            logger.warning(e)
            logger.warning(
                (
                    "Gene %s (%s:%i-%i) has no cds regions; "
                    "likely overlaps chunk border"
                ),
                gene.id,
                gene.chrom,
                gene.begin,
                gene.end,
            )
    dist = np.array(dist)
    logger.info(
        "Saw %i (%i) distances larger than double extension %i",
        np.sum(dist > dext),
        len(dist),
        dext,
    )
    out = []
    for i, g in enumerate(genes):
        if (dist[i] > dext) and (dist[i + 1] > dext):
            out.append(g)
    if len(out) == 0:
        return None
    return out


def read_genomefile(file) -> pd.DataFrame:
    """Read genome file"""
    df = read_fai_genomefile(file)
    if df is not None:
        data = df.set_index("chrom")["length"].to_dict()
        return data
    return None


@click.command(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.argument("annotation", type=click.File("r"))
@click.argument("genome", type=click.File("r"))
@click.option(
    "--window-size",
    "-w",
    type=int,
    default=1000,
    help="Window size",
)
@click.option(
    "--extension",
    "-e",
    type=int,
    default=100000,
    help=(
        "Extend windows upstream and downstream by this amount. If a"
        "neighbouring gene lies within this distance both will be discarded"
        "from analysis"
    ),
)
@click.option(
    "--output",
    "-o",
    type=click.File("w"),
    default="genomic_windows.bed",
    help="Output file",
)
@click.option(
    "gene_flag", "--gene", is_flag=True, help="Use gene instead of CDS"
)
@click.pass_context
def cli(ctx, annotation, genome, window_size, extension, output, gene_flag):  # pylint: disable=too-many-arguments,too-many-locals
    """Generate windows around genes for genic proximity analysis."""
    logger.info("Running %s %s", ctx.find_root().info_name, ctx.info_name)
    delimiter, _ = sniff_infile(annotation.name)
    compression = "gzip" if is_gzip(annotation.name) else None

    chromosomes = read_genomefile(genome)
    df = pd.read_csv(
        annotation,
        compression=compression,
        header=None,
        sep=delimiter,
        comment="#",
    )
    df.columns = [
        "chrom",
        "source",
        "feature",
        "begin",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    data = []
    for group_name, group in tqdm(df.groupby("chrom")):
        logger.info("Processing chromosome %s", group_name)
        x = process_chromosome(
            group, chromosomes[group_name], extension, gene_flag
        )
        if x is not None:
            data.extend(x)
    windows = []
    for strandedness in ["+", "-"]:
        count = 0
        logger.info("Processing %s strand", strandedness)
        for row in tqdm(data):
            if row.strand != strandedness:
                continue
            if row.strand == "+":
                upstream = np.arange(
                    row.begin,
                    row.begin - extension - window_size,
                    -window_size,
                )
                downstream = np.arange(
                    row.end, row.end + extension + window_size, window_size
                )
            else:
                upstream = np.arange(
                    row.end, row.end + extension + window_size, window_size
                )
                downstream = np.arange(
                    row.begin,
                    row.begin - extension - window_size,
                    -window_size,
                )

            count += 1
            indices = np.arange(len(upstream))
            geneid = row.id
            for begin, end, index in zip(upstream[:-1], upstream[1:], indices):
                windows.append(
                    [
                        row.chrom,
                        min(begin, end),
                        max(begin, end),
                        f"{geneid}_U{index}",
                    ]
                )
            for begin, end, index in zip(
                downstream[:-1], downstream[1:], indices
            ):
                windows.append(
                    [
                        row.chrom,
                        min(begin, end),
                        max(begin, end),
                        f"{geneid}_D{index}",
                    ]
                )

        logger.info("Processed %i genes", count)
    dfout = pd.DataFrame(windows)

    dfout.to_csv(output, sep="\t", header=False, index=False)
    logger.info("Output written to %s", output.name)
