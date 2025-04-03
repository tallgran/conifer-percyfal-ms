"""Make regions from a list of coordinates.

Make regions for features from a list of coordinates.
"""

import logging
import re
from collections import defaultdict

import click
import gffutils
import humanize
from tqdm import tqdm

from conifer.options import verbose_option

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """Make regions from a gff."""
    logger.debug("Running %s version %s", __shortname__, __version__)


@cli.command()
@click.argument("gff", type=click.Path(exists=True))
@click.argument("outfile", type=click.Path())
@verbose_option()
def create_db(gff, outfile):
    """Create database"""
    logger.info("Parsing GFF file")
    gffutils.create_db(
        gff,
        dbfn=outfile,
        force=True,
        checklines=0,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )
    logger.info("Database %s created", outfile)


@cli.command()
@click.argument("gff", type=click.Path(exists=True))
@click.argument("outfile", type=click.Path())
@click.option("--feature", "-f", help="feature to extract", default="CDS")
@click.option(
    "--annotate", "-a", help="annotate regions with feature ID", is_flag=True
)
@click.option(
    "--annotation-type",
    "-t",
    help="annotate regions with feature type",
    default="gene",
)
@click.option(
    "--min-length",
    "-l",
    help="minimum length of feature",
    default=1,
)
@verbose_option()
def make_regions(gff, outfile, feature, annotate, annotation_type, min_length):
    """Make regions of a specific feature"""
    logger.info("Parsing GFF file")
    db = gffutils.FeatureDB(gff)
    with open(outfile, "w") as fh:
        n = len(list(db.features_of_type(feature)))
        for g in tqdm(
            db.features_of_type(feature), total=n, desc="Writing regions"
        ):
            if g.end - g.start < min_length:
                logger.warning(
                    "Feature %s too short: %i-%i", g.id, g.start, g.end
                )
                continue
            if annotate:
                parents = list(db.parents(g.id, featuretype=annotation_type))
                if len(parents) > 1:
                    logger.warning("More than one parent found for %s", g.id)
                    annotation_id = parents[0].id
                elif len(parents) == 0:
                    logger.warning(
                        "No parent found for %s; using feature id", g.id
                    )
                    annotation_id = g.id
                else:
                    annotation_id = parents[0].id
                out = f"{g.chrom}\t{g.start}\t{g.end}\t{annotation_id}|{g.id}"
            else:
                out = f"{g.chrom}\t{g.start}\t{g.end}"
            print(out, file=fh)


@cli.command()
@click.argument("gffdb", type=click.Path(exists=True))
@click.argument("feature_id", type=str)
@click.option("--window-size", "-w", help="window size", default=1000)
@click.option("--feature-type", "-f", help="feature type", default="gene")
@click.option("--padding", "-p", help="padding", default=0)
@click.option(
    "--round/--no-round", "rnd", help="round to window size", default=True
)
@verbose_option()
def make_windowed_feature(
    gffdb, feature_id, window_size, feature_type, padding, rnd
):
    """Convert region to windowed regions"""
    logger.info("Parsing db file")
    db = gffutils.FeatureDB(gffdb)
    for g in db.features_of_type(feature_type):
        if g.id != feature_id:
            continue
        if rnd:
            g.start = (g.start // window_size) * window_size
            g.end = (g.end // window_size) * window_size
        g.start -= padding
        g.end += padding

        for i, pos in enumerate(range(g.start, g.end, window_size)):
            print(
                f"{g.chrom}\t{pos}\t{pos+window_size}\t{feature_id}:window{i}"
            )


@cli.command()
@click.argument("gffdb", type=click.Path(exists=True))
@click.option("--chrom", "-c", help="regex to match chromosome", default=".*")
@click.option("--feature", "-f", help="feature to extract", default=None)
@click.option(
    "--hist", "-h", help="make histogram table of lengths", is_flag=True
)
@verbose_option()
def summarize(gffdb, chrom, feature, hist):
    """Summarize GFF db file"""
    logger.info("Parsing db file")
    db = gffutils.FeatureDB(gffdb)
    features = list(db.featuretypes())
    chrom_regex = re.compile(chrom)
    logger.info("Found features %s", features)
    if feature:
        features = [feature]
        logger.info("Using feature %s", feature)
    hist_table = defaultdict(int)
    total_size = defaultdict(int)
    for feat in features:
        logger.info("Processing feature %s", feat)
        for entry in db.features_of_type(feat):
            if not chrom_regex.search(entry.chrom):
                continue
            total_size[feat] += entry.end - entry.start
            if hist:
                hist_table[entry.end - entry.start] += 1
        if hist:
            print(
                f"# Histogram of feature sizes for {feat}:"
                f"{len(hist_table)} bins"
            )
            for width, count in sorted(hist_table.items(), key=lambda x: x[0]):
                print(f"{feat}\t{width}\t{count}")
    for feat, size in total_size.items():
        print(f"Feature {feat}: {humanize.metric(size, 'bp', precision=3)}")


@cli.command()
@click.argument("gffdb", type=click.Path(exists=True))
@click.argument("feature", type=str)
@click.option(
    "--min-length", "-l", help="minimum length of feature", default=0
)
@click.option(
    "--max-length", "-m", help="maximum length of feature", default=1e6
)
@verbose_option()
def view(gffdb, feature, min_length, max_length):
    """View a feature"""
    logger.info("Parsing db file")
    db = gffutils.FeatureDB(gffdb)
    for entry in db.features_of_type(feature):
        if entry.end - entry.start < min_length:
            continue
        if entry.end - entry.start > max_length:
            continue
        print(entry)


@cli.command()
@click.argument("gffdb", type=click.Path(exists=True))
@click.option("--id", "feature_id", help="feature ID")
@click.option("--idfile", type=click.Path(exists=True), help="file with IDs")
@verbose_option()
def select(gffdb, feature_id, idfile):
    logger.info("Parsing db file %s", gffdb)
    db = gffutils.FeatureDB(gffdb)
    for entry in db.all_features():
        if feature_id:
            if entry.id == feature_id:
                print(entry)
        if idfile:
            with open(idfile) as fh:
                for line in fh:
                    if entry.id == line.strip():
                        print(entry)
