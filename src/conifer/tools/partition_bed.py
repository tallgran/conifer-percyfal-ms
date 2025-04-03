r"""Partition BEDFILE into equisized bins.

Partition BEDFILE regions into bins (subsets) of regions. Partitioning
is either done sequentially or greedily with the aim of generating
equisized partitions.

Sequential partitioning bins regions based on their index order.

Greedy partitioning first sorts regions by size and then iterates
through the remaining regions, adding to the smallest bin. Greedy
partitioning may be preferred when sequences are of very different
sizes.

"""

import logging

import click
import pybedtools

logger = logging.getLogger(__name__)


def bin_length(x):
    """Calculate total length of bin"""
    return sum(len(r) for r in x)


def cumsum(x):
    """Calculate total length of list of bins"""
    return sum(bin_length(b) for b in x)


def greedy_partition(regions, npartitions, partition):
    # Sort regions by length
    ix = sorted(
        range(len(regions)), key=lambda k: len(regions[k]), reverse=True
    )
    # Seed output regions with npartitions longest regions
    out = [[regions[i]] for i in ix[0:npartitions]]
    # Keep track of output lengths
    outlen = [len(r[0]) for r in out]
    for j in ix[npartitions : len(ix)]:  # noqa: E203
        # Get index of output set with shortest total region length
        # and add current region
        imin = outlen.index(min(outlen))
        out[imin].append(regions[j])
        outlen[imin] += len(regions[j])
    return out


def sequential_partition(regions, npartitions, partition):
    # Get total length
    length = bin_length(regions)
    binsize = length / npartitions
    out = []
    rbin = []
    size = 0
    for r in regions:
        if size >= binsize:
            out.append(rbin)
            # Adjust binsize
            if npartitions > len(out):
                binsize = (length - cumsum(out)) / (npartitions - len(out))
            rbin = []
        rbin.append(r)
        size = bin_length(rbin)
    out.append(rbin)
    return out


@click.command(help=__doc__)
@click.argument("bedfile", type=click.Path(exists=True))
@click.option(
    "partition",
    "--partition",
    "-p",
    help="requested partition number",
    default=1,
    type=int,
    show_default=True,
)
@click.option(
    "npartitions",
    "--npartitions",
    "-n",
    default=20,
    type=int,
    help="number of requested partitions",
    show_default=True,
)
@click.option(
    "--outfile", "-o", help="output file name", default=None, type=str
)
@click.option(
    "--greedy", "-g", is_flag=True, default=False, help="do greedy partition"
)
def cli(bedfile, npartitions, partition, outfile, greedy):
    regions = list(pybedtools.BedTool(bedfile))
    try:
        assert len(regions) >= npartitions
    except AssertionError:
        logger.warning(
            "Number of regions smaller than number of partitions: "
            f"'{len(regions)} < {npartitions}': "
            "lower the number of partitions with --npartitions"
        )
        raise
    try:
        assert npartitions >= partition
    except AssertionError:
        logger.warning(
            "Requested partition is larger than number of partitions: "
            f"'{partition} > {npartitions}'"
        )
        raise

    func = sequential_partition
    if greedy:
        func = greedy_partition

    out = func(regions, npartitions, partition)
    assert cumsum(out) == bin_length(regions), (
        "total bin length not equal to input bed length"
        f"{cumsum(out)} != {bin_length(regions)}"
    )
    bedout = pybedtools.BedTool(
        "\n".join(str(r) for r in out[partition - 1]), from_string=True
    )
    if outfile is not None:
        bedout.saveas(outfile)
    else:
        print(bedout)
