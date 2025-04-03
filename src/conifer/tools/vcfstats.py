"""vcfstats - calculate common population vcf stats"""

# pylint: disable=invalid-name,no-name-in-module
import itertools
import logging
import re
from pathlib import Path
from typing import Hashable

import click
import dask.array as da
import numpy as np
import sgkit as sg
import xarray as xr
from pyfastx import Fasta
from sgkit.io import vcf as sgvcf
from sgkit.model import get_contigs, num_contigs
from sgkit.stats.aggregation import (
    count_variant_alleles,
)
from sgkit.stats.popgen import diversity
from sgkit.utils import (
    conditional_merge_datasets,
    create_dataset,
    define_variable_if_absent,
)
from sgkit.variables import ArrayLikeSpec, SgkitVariables
from sgkit.window import has_windows, window_statistic
from tqdm import tqdm
from xarray import Dataset

from .. import __version__

__shortname__ = __name__.rsplit(".", maxsplit=1)[-1]

logger = logging.getLogger(__name__)
formatter = logging.Formatter(
    "%(levelname)s %(asctime)s [%(name)s:%(funcName)s]: %(message)s",
    "%Y-%m-%d %H:%M:%S",
)
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s %(asctime)s [%(name)s:%(funcName)s]: %(message)s",
)


stat_s, stat_s_spec = SgkitVariables.register_variable(
    ArrayLikeSpec(
        "stat_s",
        dims=({"windows", "variants"},),
        kind="f",
        __doc__="""Segregating sites.""",
    )
)

stat_thetaw, stat_thetaw_spec = SgkitVariables.register_variable(
    ArrayLikeSpec(
        "stat_thetaw",
        dims=({"windows", "variants"},),
        kind="f",
        __doc__="""Watterson's theta.""",
    )
)


window_coord_start, window_coord_start_spec = SgkitVariables.register_variable(
    ArrayLikeSpec(
        "window_coord_start",
        dims=("windows",),
        kind="i",
        __doc__=(
            "The physical (1-based) coordinate of the window " "start position"
        ),
    )
)


window_coord_stop, window_coord_stop_spec = SgkitVariables.register_variable(
    ArrayLikeSpec(
        "window_coord_stop",
        dims=("windows",),
        kind="i",
        __doc__=(
            "The physical (1-based, inclusive) coordinate "
            "of the window stop position"
        ),
    )
)


mask_accessible, mask_accessible_spec = SgkitVariables.register_variable(
    ArrayLikeSpec(
        "mask_accessible",
        kind="i",
        dims=({"contigs", "windows"},),
        __doc__="""Accessible sites in each contig / window.""",
    )
)


mask_inaccessible, mask_inaccessible_spec = SgkitVariables.register_variable(
    ArrayLikeSpec(
        "mask_inaccessible",
        kind="i",
        dims=({"contigs", "windows"},),
        __doc__="""Inccessible sites in each contig / window.""",
    )
)


def validate_sample_set(ctx, param, value):  # pylint: disable=unused-argument
    """Parse sample set parameter"""
    retval = []
    for item in value:
        if re.search(r":", item) is not None:
            m = re.match(
                r"^(?P<label>[^:,\s]*):?(?P<ss>[^:,\s]+),?(?P<mask>[^:,\s]*)$",
                item,
            )
        else:
            m = re.match(r"^(?P<ss>[^:,\s]+),?(?P<mask>[^:,\s]*)$", item)
        try:
            retval.append(tuple(m.groupdict().values()))
        except AttributeError as exc:
            raise click.BadParameter(
                f"bad format '{item}'\n"
                "format must be one of \n"
                "  <sample_set_file> \n"
                "  <sample_set_file>,<mask_file>\n"
                "  <sample_set_label>:<sample_set_file>\n"
                "  <sample_set_label>:<sample_set_file>,<mask_file>"
            ) from exc
    return retval


def mask2bool(x: str):
    """Convert pyfastx mask sequence to numpy boolean array"""
    y = np.frombuffer(x.encode("ascii"), dtype="u1") - 48
    return np.array(y, dtype=bool)


def add_window_coord(  # pylint: disable=too-many-locals
    ds: Dataset,
    *,
    window_size: int,
    merge: bool = True,
) -> Dataset:
    """Add start and stop coordinates of windows translated from
    variant indices to Dataset. Modeled on
    sgkit.window._window_per_contig."""
    n_contigs = num_contigs(ds)
    contig_ids = np.arange(n_contigs)
    contig_lengths = ds.get("contig_length", [None] * n_contigs)
    window_contig = ds["window_contig"]
    window_stop_index = (
        np.searchsorted(window_contig, contig_ids, side="right") - 1
    )
    window_upper_bound = ds["variant_position"][
        ds["window_stop"][window_stop_index] - 1
    ]
    n_windows_by_contig = np.bincount(ds["window_contig"].values)

    contig_window_coord_contigs = []
    contig_window_coord_starts = []
    contig_window_coord_stops = []

    for i, _ in enumerate(get_contigs(ds)):
        length = contig_lengths[i]
        nwin = n_windows_by_contig[i]
        if length is None:
            length = window_upper_bound[i] + window_size
        windows = np.arange(1, length, window_size)
        starts = windows[:nwin]
        stops = windows[1 : (nwin + 1)] - 1
        if len(stops) < len(starts):
            stops = np.concatenate((stops, [length]))
        contig_window_coord_starts.append(starts)
        contig_window_coord_stops.append(stops)
        contig_window_coord_contigs.append(np.full_like(starts, i))

    window_coord_starts = np.concatenate(contig_window_coord_starts)
    window_coord_stops = np.concatenate(contig_window_coord_stops)
    new_ds = create_dataset(
        {
            window_coord_start: (
                "windows",
                window_coord_starts,
            ),
            window_coord_stop: (
                "windows",
                window_coord_stops,
            ),
        }
    )
    return conditional_merge_datasets(ds, new_ds, merge)


def add_mask_accessibility(  # pylint: disable=too-many-locals
    ds: Dataset,
    *,
    mask_faidx: Fasta = None,
    window_size: int = None,
    merge: bool = True,
):
    """Add mask summary for normalizing statistics.

    Note: non-window stastics are calculated per site whereby mask summaries
    are reported per contig.
    """
    if mask_faidx is None:
        return ds
    logger.info("add mask accessibility to dataset")
    n_contigs = num_contigs(ds)
    if "windows" in ds.sizes:
        n_windows_by_contig = np.bincount(ds["window_contig"].values)
    else:
        n_windows_by_contig = np.repeat(1, n_contigs)

    contig_accessible = []
    contig_inaccessible = []

    for i, contig in enumerate(get_contigs(ds)):
        logger.info("add mask accessibility for %i, %s", i, contig)
        x = mask2bool(mask_faidx[str(contig)].seq)
        contig_length = len(x)
        if window_size is None:
            macc = (~x,)
            minacc = (x,)
        else:
            window_start = np.arange(window_size, contig_length, window_size)
            macc = np.array_split(~x, window_start)[: n_windows_by_contig[i]]
            minacc = np.array_split(x, window_start)[: n_windows_by_contig[i]]
        contig_accessible.append([np.sum(m) for m in macc])
        contig_inaccessible.append([np.sum(m) for m in minacc])
    mask_accessible_sites = np.concatenate(contig_accessible)
    mask_inaccessible_sites = np.concatenate(contig_inaccessible)
    if window_size is None:
        mask_dim = "contigs"
    else:
        mask_dim = "windows"
    new_ds = create_dataset(
        {
            mask_accessible: (
                mask_dim,
                mask_accessible_sites,
            ),
            mask_inaccessible: (
                mask_dim,
                mask_inaccessible_sites,
            ),
        }
    )
    logger.info(
        "Added mask information for %i (%i)",
        len(mask_accessible_sites),
        mask_dim,
    )
    return conditional_merge_datasets(ds, new_ds, merge)


def theta_s(
    ds: Dataset,
    *,
    variant_allele_count: Hashable = sg.variables.variant_allele_count,
    stat_diversity: Hashable = sg.variables.stat_diversity,
    merge: bool = True,
) -> Dataset:
    """Compute Watterson's theta and number of segregating sites.

    Based on sgkit.popgen.Tajimas_D, which calculates these statistics
    internally without saving them. Note that this calculation does
    *not* work for cohorts!
    """
    ds = define_variable_if_absent(
        ds,
        sg.variables.variant_allele_count,
        variant_allele_count,
        count_variant_alleles,
    )
    ds = define_variable_if_absent(
        ds, sg.variables.stat_diversity, stat_diversity, diversity
    )
    sg.variables.validate(
        ds,
        {
            variant_allele_count: sg.variables.variant_allele_count_spec,
            stat_diversity: sg.variables.stat_diversity_spec,
        },
    )

    ac = ds[variant_allele_count]
    ac = da.asarray(ac)

    # count segregating. Note that this uses the definition in tskit,
    # which is the number of alleles - 1. In the biallelic case this
    # gives us the number of non-monomorphic sites.
    s = (ac > 0).sum(axis=1) - 1
    # Segregating sites less than 0 make no sense
    s = np.array(s, dtype="float")
    s[np.where(s < 0)] = np.nan

    if has_windows(ds):
        s = window_statistic(
            s,
            np.sum,
            ds.window_start.values,
            ds.window_stop.values,
            dtype=s.dtype,
            axis=0,
        )

    # assume number of chromosomes sampled is constant for all variants
    # NOTE: even tho ac has dtype uint, we promote the sum to float
    #       because the computation below requires floats
    n = ac.sum(axis=1, dtype="float").max()

    # (n-1)th harmonic number
    a1 = (1 / da.arange(1, n)).sum()

    # calculate Watterson's theta (absolute value)
    theta = s / a1
    vname = "variants"
    if has_windows(ds):
        vname = "windows"
    new_ds = create_dataset(
        {
            stat_s: ((vname,), s),
            stat_thetaw: ((vname,), theta),
        }
    )
    return conditional_merge_datasets(ds, new_ds, merge)


def load_mask(mask):
    """Load mask file."""
    if mask is None:
        return None
    logger.info("loading mask file %s", mask)
    return Fasta(mask)


def load_samples(*, sampleset, sample_ids):
    """Load sample set names and mask files.

    Return tuple of sampleset label, list of samples, and mask file.
    """
    value = []
    if len(sampleset) == 0:
        return (("ALL", sample_ids, None),)
    for name, ss, mask in sampleset:
        if name is None:
            name = Path(ss).name
        samples = [
            x for x in Path(ss).read_text(encoding="utf-8").split("\n") if x
        ]
        mask_faidx = load_mask(mask)
        value.append((name, samples, mask_faidx))
    return tuple(
        value,
    )


def calculate_stats(ds, stats):  # pylint: disable=redefined-outer-name
    """Calculate statistics."""
    if "pi" in stats:
        logger.info("Calculating pi...")
        ds = sg.diversity(ds, merge=True)
    if "TajD" in stats:
        logger.info("Calculating Tajimas D...")
        ds = sg.Tajimas_D(ds, merge=True)
    if ("S" in stats) or ("thetaW" in stats):
        logger.info("Calculating S and Watterson's theta...")
        ds = theta_s(ds, merge=True)
    if "Fst" in stats:
        logger.info("Calculating Fst...")
        ds = sg.Fst(ds, merge=True)
    if "divergence" in stats:
        logger.info("Calculating dxy...")
        ds = sg.divergence(ds, merge=True)
    return ds


def subset_by_regions(ds):
    """Subset by regions. NB: currently only subsets by chromosome.
    Needed to prevent pathological case where there are zero variants
    on a chromosome which causes windowing to fail."""
    logger.info("Subsetting dataset to contigs with variants")
    contigs = list(set(ds.variant_contig.values))
    return ds.isel(contigs=contigs)


def subset_by_mask(ds: Dataset, *, mask_faidx: Fasta = None, offset: int = 1):
    """Subset dataset by mask. The mask is accessed through the pyfaidx API"""
    if mask_faidx is None:
        return ds
    logger.info("Subsetting dataset by mask in %s", mask_faidx.file_name)
    contig_count = np.bincount(
        ds.variant_contig.values, minlength=len(get_contigs(ds))
    )
    variant_index = []
    for i, contig in tqdm(enumerate(get_contigs(ds))):
        logger.debug("Subsetting mask on contig %s", contig)
        if contig_count[i] == 0:
            continue
        index = ds.variants.values[ds.variant_contig == i]
        x = mask2bool(mask_faidx[str(contig)].seq)
        variant_pos = ds.variant_position[index] - offset
        j = x[variant_pos]
        variant_index.append(index[~j])
    i = np.concatenate(variant_index)
    logger.info(
        "selected %.3f (%.3f%%) variants",
        len(i) / len(ds.variants),
        len(i) / len(ds.variants) * 100,
    )
    return ds.sel(variants=i)


def set_multiway_cohort(
    ds: Dataset,
    ss1: list,  # pylint: disable=unused-argument
    ss2: list,
) -> Dataset:
    """Set multiway cohort."""
    sample_cohort = xr.DataArray(
        np.full(len(ds.samples), 0, dtype="i4"), dims="samples"
    )
    i2 = ds.sample_id.isin(ss2).values
    sample_cohort[i2] = 1
    ds["sample_cohort"] = sample_cohort
    return ds


def load_dataset(path):
    """Load dataset from path."""
    if path is None:
        return None
    logger.info("Loading path %s", path)
    ds = sg.load_dataset(path)
    # Fallback on defaults
    if "sample_cohort" not in ds.keys():
        sample_cohort = xr.DataArray(
            np.full(len(ds.samples), 0, dtype="i4"), dims="samples"
        )
        ds["sample_cohort"] = xr.DataArray(sample_cohort, dims="samples")
    return ds


def remove_monomorphic_sites(
    ds: Dataset,
    *,
    variant_allele_count: Hashable = sg.variables.variant_allele_count,
):
    """Remove monomorphic sites from dataset."""
    logger.info("Removing monomorphic sites...")
    ac = sg.count_variant_alleles(ds)[variant_allele_count]
    s = (ac > 0).sum(axis=1) - 1
    variants_sel = np.where(s > 0)[0]
    logger.info(
        "Removed  %.3f (%.3f%%) monomorphic / missing data sites",
        (len(s) - len(variants_sel)) / len(s),
        (len(s) - len(variants_sel)) / len(s) * 100,
    )
    return ds.sel(variants=variants_sel)


@click.group(help=__doc__, name=__shortname__)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    """vcfstats - calculate common population vcf stats"""
    ctx.ensure_object(dict)
    ctx.obj["VERSION"] = __version__


@cli.command()
@click.argument("zarrpath")
@click.argument("vcf", type=click.Path(exists=True), nargs=-1)
@click.option(
    "--regions", type=str, default=None, required=False, multiple=True
)
@click.option("--tempdir", default="/tmp", type=click.Path())
@click.pass_context
def vcf_to_zarr(ctx, zarrpath, vcf, regions, tempdir):
    """Convert vcf to zarr"""
    tempdir = Path(tempdir) / Path(zarrpath).name
    logger.info("Running %s vcf-to-zarr", ctx.find_root().info_name)
    if len(regions) > 0:
        regions = [r.split(",") for r in regions]
    else:
        regions = None
    logger.info("Converting vcf files %s to zarr: %s", ",".join(vcf), zarrpath)
    sgvcf.vcf_to_zarr(vcf, zarrpath, regions=regions, tempdir=str(tempdir))
    logger.info("Conversion to zarr done")


@cli.command()
@click.argument("path", type=click.Path(exists=True))
@click.argument("regions", type=str, default=None, required=False, nargs=-1)
@click.option(
    "--outfile-prefix",
    "-o",
    default="vcfstats",
    type=str,
    help="outfile prefix",
)
@click.option(
    "-s",
    "--stats",
    default=["pi", "thetaW", "S", "TajD"],
    type=click.Choice(["pi", "thetaW", "S", "TajD", "Fst", "dxy"]),
    multiple=True,
    help=(
        "Apply oneway and/or multiway statistics to data "
        "Oneway statistics (pi, thetaW, S, and TajD) are "
        "calculated separately for every sampleset. Multiway "
        "statistics (Fst, dxy) are calculated for all pairwise "
        "sampleset combinations."
    ),
)
@click.option(
    "--sample-set",
    default=None,
    type=click.UNPROCESSED,
    callback=validate_sample_set,
    multiple=True,
    required=False,
    help=(
        "The sample set specification. The sample set consists of a file with "
        "sample ids, one per line, without header. By default, the file "
        "basename will be used as sample set label. Optionally, a label "
        "can be assigned to to the "
        "sample set as <sample_set_label>:<sample_set_file>. In addition, "
        "a sample set specific mask can be specified by adding a comma "
        "separator and file name as "
        "<sample_set_label>:<sample_set_file>,<mask_file>. "
        "The mask file is a 0/1-encoded fasta file, where 1=discard, 0=keep. "
        "The option can be applied multiple times. For two way statistics all "
        "pairwise comparisons between sample sets will be made. "
    ),
)
@click.option(
    "-w",
    "--window-size",
    default=None,
    type=int,
    help=("Window size"),
)
@click.pass_context
def stats(  # pylint: disable=too-many-arguments,too-many-locals
    ctx,
    path,
    regions,
    outfile_prefix,
    stats,  #  pylint: disable=redefined-outer-name
    sample_set,
    window_size,
    **kwargs,
):  # noqa A001, pylint: disable=unused-argument,too-many-arguments,too-many-branches,too-many-statements
    """Calculate statistics on variant data."""
    logger.info("Running %s stats", ctx.find_root().info_name)
    oneway = {
        "pi": "stat_diversity",
        "thetaW": "stat_thetaW",
        "S": "stat_S",
        "TajD": "stat_Tajimas_D",
    }
    multiway = {
        "Fst": "stat_Fst",
        "dxy": "stat_divergence",
    }
    ds = load_dataset(path)
    if window_size is None:
        cols = ["contig_id", "variant_position", "variant_contig"]
    else:
        cols = [
            "contig_id",
            "window_coord_start",
            "window_coord_stop",
            "window_contig",
        ]
    sample_set = load_samples(
        sampleset=sample_set,
        sample_ids=ds.sample_id.values.tolist(),
    )
    # Oneway statistics
    _stats = list(set(stats).intersection(oneway.keys()))
    if len(_stats) == 0:
        logger.info("No oneway statistics defined; skipping")
    else:
        for name, ss, mask_faidx in sample_set:
            logger.info("Processing sample set %s; %i samples", name, len(ss))
            outcols = cols
            x = ds.sel(samples=ds.sample_id.isin(ss).values)
            x = subset_by_mask(x, mask_faidx=mask_faidx)
            x = remove_monomorphic_sites(x)
            x = subset_by_regions(x)
            if window_size is not None:
                logger.info("adding windows by position, size %i", window_size)
                x = sg.window_by_position(x, size=window_size, offset=1)
                x = add_window_coord(x, window_size=window_size)
            x = add_mask_accessibility(
                x, mask_faidx=mask_faidx, window_size=window_size
            )
            x = calculate_stats(x, _stats)
            outprefix = f"{outfile_prefix}-{name}-oneway"
            if window_size is not None:
                outprefix = f"{outprefix}-w{window_size}"
            if mask_accessible in x:
                outcols = outcols + [mask_accessible, mask_inaccessible]
            store = f"{outprefix}.zarr"
            outcols = outcols + list(map(oneway.get, _stats))
            try:
                logger.info("Saving results to zarr store %s", store)
                sg.save_dataset(
                    x[outcols].unify_chunks(), store, auto_rechunk=True
                )
            except Exception as e:
                logger.error("error: %s", e)
                raise

    # Multi-way statistics
    _stats = list(set(stats).intersection(multiway.keys()))
    outcols = cols + list(map(multiway.get, stats))
    if (len(sample_set) < 2) or (len(_stats) == 0):
        logger.info(
            "No multiway statistics defined and/or only one sampleset; "
            "skipping multiway statistics."
        )
    else:
        outcols = cols
        for ss1, ss2 in itertools.combinations(sample_set, 2):
            name1, samples1, mask_faidx1 = ss1
            name2, samples2, mask_faidx2 = ss2
            samples = samples1 + samples2
            logger.info(
                "Processing sample sets %s and %s; %i samples",
                name1,
                name2,
                len(samples),
            )
            x = ds.sel(samples=ds.sample_id.isin(samples).values)
            x = subset_by_mask(x, mask_faidx=mask_faidx1)
            x = subset_by_mask(x, mask_faidx=mask_faidx2)
            x = remove_monomorphic_sites(x)
            x = subset_by_regions(x)
            x = set_multiway_cohort(x, samples1, samples2)
            if window_size is not None:
                logger.info("adding windows by position, size %i", window_size)
                x = sg.window_by_position(x, size=window_size, offset=1)
                x = add_window_coord(x, window_size=window_size)
            x = add_mask_accessibility(
                x, mask_faidx=mask_faidx, window_size=window_size
            )
            x = calculate_stats(x, _stats)
            outprefix = f"{outfile_prefix}-{name1}-{name2}-multiway"
        if window_size is not None:
            outprefix = f"{outprefix}-w{window_size}"
        if mask_accessible in x:
            outcols = outcols + [mask_accessible, mask_inaccessible]
        store = f"{outprefix}.zarr"
        outcols = outcols + list(map(multiway.get, _stats))
        try:
            logger.info("Saving results to zarr store %s", store)
            sg.save_dataset(
                x[outcols].unify_chunks(), store, auto_rechunk=True
            )
        except Exception as e:
            logger.error("error: %s", e)
            raise


@click.command()
@click.argument("zarrpath", type=click.Path(exists=True))
@click.argument(
    "statistic",
    type=click.Choice(["stat_diversity", "stat_TajimasD"]),
    nargs=-1,
)
@click.option(
    "--function",
    "-f",
    default=["np.sum", "np.mean", "np.var"],
    type=click.Choice(["np.sum", "np.var", "np.mean", "np.median"]),
    multiple=True,
    help=("Apply function to numerical quantity in zarr path"),
)
@click.pass_context
def summary(ctx, zarrpath, statistic):  # pylint: disable=unused-argument
    """Summarize statistics in zarrpath"""
