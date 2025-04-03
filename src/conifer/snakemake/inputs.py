"""Snakemake input functions"""

import glob
import pathlib
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from conifer.snakemake import samples
from snakemake.io import expand

BIALLELIC_VCF = glob.glob(
    "data/vcf/PA_chr[0-9][0-9]_[0-9][0-9].vcf.gz"
) + glob.glob("data/vcf/PA_chr[0-9][0-9]_[0-9].vcf.gz")
BIALLELIC_VCF = [
    pathlib.Path(x).name.rstrip(r"\.vcf.gz") for x in BIALLELIC_VCF
]
CHROMOSOMES = [f"PA_chr{x:02d}" for x in range(1, 13)]
CHROMOSOME_MAP = defaultdict(list)
for pfx in BIALLELIC_VCF:
    KEY = "_".join(pfx.split("_")[:2])
    CHROMOSOME_MAP[KEY].append(pfx)


def bam_files(wildcards, samplelist):
    """Retrieve list of bam files"""

    def _format(x):
        return f"{wildcards.data}/bam/dedup.ir.{x}.bam"

    begin = 0
    end = len(samplelist)
    if "samplepartition" in wildcards.keys():
        begin = int((wildcards.samplepartition) - 1) * int(
            wildcards.sample_size
        )
    if "sample_size" in wildcards.keys():
        end = begin + int(wildcards.sample_size)

    return samplelist.loc[begin:end].SM.map(_format).tolist()


def bam_files_population_sampleset(wildcards):
    """Retrieve list of bam files for wildcards.sampleset for a given
    population type"""

    def _format(x):
        return f"data/bam/dedup.ir.{x}.bam"

    sample_df = samples.read_sampleset(
        wildcards,
        fmt="samples-{populationtype}-{sampleset}-of-{nsamplesets}.tsv",
    )
    return sample_df.SM.map(_format).tolist()


def d4_from_sampleset_input(wildcards):
    """Return d4 files for a sampleset"""
    path = f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.MQ}"

    def _format(x):
        return f"{path}/{x}.per-base.d4"

    sample_df = samples.read_sampleset(wildcards)
    return sample_df.SM.map(_format).tolist()


def d4_count_min_coverage_chunk_input(wildcards):
    """Retrieve d4 files to run count min coverage in chunks"""
    path = f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.MQ}"

    def _format(x):
        return f"{path}/{x}.per-base.d4"

    sample_df = samples.read_sampleset(wildcards)
    return sample_df.SM.map(_format).tolist()


def d4_count_min_coverage_aggregate_input(wildcards):
    """Retrieve d4 chunk files to merge for wildcards.sampleset"""
    path = (
        f"{wildcards.data}/aggregate_coverage/"
        f"{wildcards.sampleset}/"
        f"MQ{wildcards.MQ}/chunks/"
        f"count_ge{wildcards.coverage}."
        f"{{partition}}-of-{wildcards.npartitions}.d4"
    )
    return [path.format(partition=x) for x in list(range(1, 101, 1))]


def d4_sum_coverage_aggregate_input(wildcards):
    """Retrieve d4 chunk files to merge for wildcards.sampleset"""
    path = (
        f"{wildcards.data}/aggregate_coverage/"
        f"{wildcards.sampleset}/"
        f"MQ{wildcards.MQ}/chunks/"
        f"sum.{{partition}}-of-{wildcards.npartitions}.d4"
    )
    return [path.format(partition=x) for x in list(range(1, 101, 1))]


def partition_samples_on_coverage_input(wildcards):
    """Return mosdepth summary files for partitioning samples on coverage"""
    sample_df = pd.read_csv(wildcards.samplesheet, sep="\t")
    data = re.sub("resources", "data", wildcards.resources)

    def _format(x):
        fmt = f"{data}/mosdepth_coverage/MQ10/{x}.mosdepth.summary.txt"
        return fmt

    return sample_df.SM.map(_format).tolist()


def bedtools_unionbedgraph_input(wildcards, ext="", samplenames=False):
    """Retrieve bedgraph files for a sampleset"""
    path = f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.mapping_quality}"

    def _format(x):
        fmt = (
            f"{path}/{wildcards.partition_method}/{wildcards.npartitions}/"
            f"{x}.{wildcards.partition}-of-{wildcards.npartitions}."
            f"per-base.bed.gz{ext}"
        )
        return fmt

    sample_df = samples.read_sampleset(wildcards)

    if samplenames:
        return sample_df.SM.tolist() if len(sample_df) > 0 else []

    return sample_df.SM.map(_format).tolist()


def d4_concat_input(wildcards):
    """Retrieve d4 files for wildcards.sampleset"""
    fmt = (
        f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.mapping_quality}/"
        f"{wildcards.sampleset}/{wildcards.partition_method}/"
        f"{{i}}-of-{wildcards.npartitions}.per-base.sum.d4"
    )
    ilist = np.arange(int(wildcards.npartitions)) + 1
    return [fmt.format(i=i) for i in ilist]


def get_resources_path(wildcards):
    """Retrieve resources path"""
    return pathlib.Path("resources") / pathlib.Path(
        wildcards.data
    ).relative_to("data")


def samtools_depth_merge_subset_samples_input(
    wildcards,
):  # pylint: disable=unused-argument
    """Retrieve samtools depth files for a sampleset"""


def bedgraph_concat_min_individual_cutoff_inputs(
    wildcards,
):  # pylint: disable=unused-argument
    """Retrieve bedgraph files for a sampleset"""


def _vcf_to_zarr(chrom, ext=""):
    """Retrieve vcf or tbi files for a chromosome"""
    fmt = f"data/vcf/{{chrom}}.biallelic.vcf.gz{ext}"
    return [fmt.format(chrom=x) for x in CHROMOSOME_MAP[chrom]]


def vcf_to_zarr_vcf_input(wildcards):
    """Retrieve vcf files for a chromosome"""
    return _vcf_to_zarr(wildcards.chromosome)


def vcf_to_zarr_tbi_input(wildcards):
    """Retrieve tbi files for a chromosome"""
    return _vcf_to_zarr(wildcards.chromosome, ext=".tbi")


def combine_vcfstats_oneway_statistics(wildcards):
    """Retrieve vcfstats files for all samples"""
    fmt = (
        f"{wildcards.results}/diversity.vcfstats/{wildcards.filter}/"
        f"{wildcards.sampleset}/{wildcards.chromosome}/{{pfx}}.{wildcards.source}."
        f"{wildcards.feature}{wildcards.covmask}-{wildcards.sampleset}-"
        f"oneway{wildcards.window}.zarr"
    )
    return [fmt.format(pfx=x) for x in CHROMOSOME_MAP[wildcards.chromosome]]


def summarize_vcftools_statistics_masks(wildcards, config, tbi=False):
    """Make vcftools statistics summary input for a given mask"""
    results = []
    is_sample = config["mask"]["sampleset"][wildcards.sampleset].get(
        "is_sample", False
    )
    label = "" if is_sample else ".p100"
    fmt = (
        f"data/mask/{wildcards.sampleset}/"
        f"{{label}}.filter_min{{min}}_max{{max}}{label}.bed.gz"
    )
    if tbi:
        fmt += ".tbi"
    masks = config["mask"]["sampleset"][wildcards.sampleset]
    for key, params in masks.items():
        if key == "window_size" or key == "is_sample":
            continue
        results.append(fmt.format(**params))
    return results


def summarize_vcftools_statistics_features(wildcards, config, tbi=False):  # pylint: disable=unused-argument
    """Make vcftools statistics summary input for a given feature"""
    results = []
    fmt = "data/resources/regions/{feature}.bed.gz"
    if tbi:
        fmt += ".tbi"
    features = config["mask"]["features"]
    for feature in features:
        results.append(fmt.format(feature=feature))
    return results


def multiqc_coverage_input(wildcards, samplelist):
    """Retreive multiqc coverage inputs"""
    if not isinstance(samplelist, pd.Series):
        samplelist = pd.Series(samplelist)

    def _stat_format(x, stat="hist"):
        fmt = (
            f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.mapping_quality}/"
            f"{x}.per-base.d4.{stat}"
        )
        return fmt

    def _mosdepth_format(x, stat="global.dist"):
        fmt = (
            f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.mapping_quality}/"
            f"{x}.mosdepth.{stat}.txt"
        )
        return fmt

    def _samtools_stat_format(x):
        fmt = f"results/samtools_stats/{x}.txt"
        return fmt

    ret = (
        samplelist.apply(_stat_format).tolist()
        + samplelist.apply(_stat_format, args=("mean",)).tolist()
        + samplelist.apply(_mosdepth_format).tolist()
        + samplelist.apply(_mosdepth_format, args=("summary",)).tolist()
        + samplelist.apply(_samtools_stat_format).tolist()
    )
    return ret


def multiqc_qc_input(config):
    """Retrieve multiqc files for all qc analyses."""
    inputs = expand(
        "data/mosdepth_coverage/MQ{MQ}/multiqc_input.txt",
        MQ=config["coverage"]["MQ"],
    )

    return inputs


def multiqc_configfile(wildcards, output):  # pylint: disable=unused-argument
    """Retrieve multiqc config file"""
    configfile = pathlib.Path(output.html.replace(".html", ".yaml"))
    if configfile.exists():
        return configfile
    return ""
