"""Presence / absence analysis workflow.

Compile presence / absence data based on coverage data. The output is
a matrix of regions (rows) and samples (columns) encoded by their
coverage or as 0/1 for absence/precsence. Regions could be features
(e.g., genes) or single base pair positions with 1 indicating presence
and 0 indicating absence. Matrices are used as inputs for making PCA
plots.

"""
import re
import sys
from conifer.snakemake.config import SchemaFiles
from pathlib import Path
from snakemake.utils import validate
import pandas as pd

from conifer.snakemake.lmod import get_envmodules


envvars:
    "SCRATCH",


PICAB_TE_ANNOTATION = "Picab02_230926_at01_longest"
PICAB_ANNOTATION = "Picab02_230926_at01_longest_no_TE_sorted"
PICAB_PSEUDOGENE_ANNOTATION = "Picab02_230926_at01_longest_no_TE_sorted.pseudogene"
SCRATCH = Path(os.path.normpath(os.environ["SCRATCH"]))


##############################
# Configuration
##############################
configfile: Path("config/config.yaml")


validate(config, schema=SchemaFiles.CONFIGURATION_SCHEMA)

samples = pd.read_table(config["samples"])
validate(samples, schema=SchemaFiles.SAMPLES_SCHEMA)


wildcard_constraints:
    data="data",
    interim="data/interim",
    reports="reports",
    sample=f"({'|'.join(samples.SM.tolist())})",


ALL = {
    "presabs": expand(
        "data/presabs/{feature}/{sample}.sum.tsv",
        feature=["gene", "CDS", "intron", "pseudogene", "TE", "random"],
        sample=samples.SM.tolist(),
    )
}


rule all:
    input:
        **ALL,


rule all_d4:
    input:
        expand(
            SCRATCH / "data/mosdepth_coverage/MQ10/{sample}.per-base.d4",
            sample=samples.SM.tolist(),
        ),


rule all_d4_link:
    input:
        expand(
            "data/mosdepth_coverage/MQ10/{sample}.per-base.d4",
            sample=samples.SM.tolist(),
        ),


##############################
# Atomic rules
##############################
rule gffutils_make_regions:
    """Make regions from GFF database."""
    output:
        bed="data/resources/presabs/{feature}.bed",
    input:
        db=f"data/resources/{PICAB_ANNOTATION}.liftover.gff3.db",
    conda:
        "../envs/bedtools.yaml"
    envmodules:
        *config["envmodules"]["bedtools"],
    wildcard_constraints:
        feature="gene|CDS|exon|intron|UTR",
    benchmark:
        "benchmarks/gffutils_make_regions/data/resources/presabs/{feature}.bed.benchmark.txt"
    log:
        "logs/gffutils_make_regions/data/resources/presabs/{feature}.bed.log",
    threads: 1
    shell:
        """conifer-gffutils make-regions {input.db} {output.bed} -a --feature {wildcards.feature}"""


rule gffutils_make_TE_regions:
    """Special case rule for TE annotations."""
    output:
        bed="data/resources/presabs/TE.bed",
    input:
        db=f"data/resources/{PICAB_TE_ANNOTATION}.TE.liftover.gff3.db",
    envmodules:
        *config["envmodules"]["bedtools"],
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/gffutils_make_TE_regions/data/resources/presabs/TE.bed.benchmark.txt"
    log:
        "logs/gffutils_make_TE_regions/data/resources/presabs/TE.bed.log",
    threads: 1
    shell:
        """conifer-gffutils make-regions {input.db} {output.bed} -a --feature gene"""


rule gffutils_make_pseudogene_regions:
    """Special case rule for pseudogene annotations."""
    output:
        bed="data/resources/presabs/pseudogene.bed",
    input:
        db=f"data/resources/{PICAB_PSEUDOGENE_ANNOTATION}.liftover.gff3.db",
    envmodules:
        *config["envmodules"]["bedtools"],
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/gffutils_make_pseudogene_regions/data/resources/presabs/pseudogene.bed.benchmark.txt"
    log:
        "logs/gffutils_make_pseudogene_regions/data/resources/presabs/pseudogene.bed.log",
    threads: 1
    shell:
        """conifer-gffutils make-regions {input.db} {output.bed} -a --feature gene"""


rule make_pabies_chromosome_file:
    """Make sorting file for bedtools sort."""
    output:
        "data/resources/pabies-2.0.fa.chromosomes.txt",
    input:
        "data/resources/pabies-2.0.fa.fai",
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/make_pabies_chromosome_file/data/resources/pabies-2.0.fa.chromosomes.txt.benchmark.txt"
    log:
        "logs/make_pabies_chromosome_file/data/resources/pabies-2.0.fa.chromosomes.txt.log",
    threads: 1
    shell:
        """cut -f 1 {input} > {output}"""


rule shuffle_pseudogene_regions:
    """Shuffle pseudogene regions by modifying the chromosome names."""
    output:
        bed="data/resources/presabs/random.bed",
    input:
        bed="data/resources/presabs/pseudogene.bed",
        txt="data/resources/pabies-2.0.fa.chromosomes.txt",
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/shuffle_pseudogene_regions/data/resources/presabs/random.bed.benchmark.txt"
    log:
        "logs/shuffle_pseudogene_regions/data/resources/presabs/random.bed.log",
    threads: 1
    shell:
        """
        paste <(cut -f 1 {input.bed}) <(cut -f 2,3 {input.bed} | shuf --random-source=<(yes 42)) |\
        bedtools sort -g {input.txt} |\
        paste - <(cut -f 4 {input.bed}) > {output.bed}
        """


rule mosdepth_d4:
    output:
        perbase=str(
            SCRATCH / "data/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4"
        ),
    input:
        bam="data/bam/dedup.ir.{sample}.bam",
        bai="data/bam/dedup.ir.{sample}.bam.bai",
    params:
        prefix=lambda wildcards: SCRATCH
        / f"data/mosdepth_coverage/MQ{wildcards.mapping_quality}/{wildcards.sample}",
    threads: 10
    resources:
        mem_mb=60000,
    conda:
        "../envs/mosdepth.yaml"
    envmodules:
        *get_envmodules(config, "mosdepth"),
    log:
        "logs/data/mosdepth_coverage/MQ{mapping_quality}/{sample}.mosdepth.log",
    priority: 100
    benchmark:
        "logs/benchmark/data/mosdepth_coverage/MQ{mapping_quality}/{sample}.mosdepth.benchmark.txt"
    shell:
        "mosdepth --d4 -x {params.prefix} -Q {wildcards.mapping_quality} {input.bam} -t {threads}"


rule link_d4:
    """Link d4 from scratch to data directory."""
    output:
        d4="data/mosdepth_coverage/MQ10/{sample}.per-base.d4",
    input:
        d4=lambda wildcards: storage(
            "file://"
            + str(
                SCRATCH / f"data/mosdepth_coverage/MQ10/{wildcards.sample}.per-base.d4"
            )
        ),
    conda:
        "../envs/storage.yaml"
    benchmark:
        "benchmarks/link_d4/data/mosdepth_coverage/MQ10/{sample}.per-base.d4.benchmark.txt"
    log:
        "logs/link_d4/data/mosdepth_coverage/MQ10/{sample}.per-base.d4.log",
    threads: 1
    shell:
        """cp -d {input.d4} {output.d4}"""


rule get_individual_coverage_by_feature:
    """Get the individual coverage based on a feature"""
    output:
        bed=temp("data/presabs/{feature}/{sample}.bed"),
    input:
        bed="data/resources/presabs/{feature}.bed",
        d4="data/mosdepth_coverage/MQ10/{sample}.per-base.d4",
    envmodules:
        *config["envmodules"]["bedtools"],
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/get_individual_coverage/data/presabs/{feature}/{sample}.bed.benchmark.txt"
    log:
        "logs/get_individual_coverage/data/presabs/{feature}/{sample}.bed.log",
    threads: 1
    shell:
        """d4tools view -R {input.bed} {input.d4} |\
        bedtools intersect -a - -b {input.bed} -wb |\
        csvtk cut -t -T -f 1,2,3,8,4 > {output.bed}"""


rule bgzip_compress_and_index:
    """Run bgzip compress and index"""
    output:
        gz="data/presabs/{feature}/{sample}.bed.gz",
        gzi="data/presabs/{feature}/{sample}.bed.gz.gzi",
    input:
        gz="data/presabs/{feature}/{sample}.bed",
    envmodules:
        *config["envmodules"]["tabix_bgzip"],
    conda:
        "../envs/tabixbgzip.yaml"
    benchmark:
        "benchmarks/bgzip_index/data/presabs/{feature}/{sample}.bed.gz.gzi.benchmark.txt"
    log:
        "logs/bgzip_index/data/presabs/{feature}/{sample}.bed.gz.gzi.log",
    threads: 1
    shell:
        """bgzip -i {input.gz}"""


## Somehow sum coverage over gene; the normalization will be taken
## care of later by post-processing script. Coverage below a certain
## threshold is set to zero, above to 1. The output should therefore
## be a weighted average over the feature, representing proportion of
## missingness over the feature body
rule summarize_coverage_over_feature:
    """Summarize coverage over feature.

    Recode input to 0/1 and compile the weighted output sum.

    """
    output:
        tsv="data/presabs/{feature}/{sample}.sum.tsv",
    input:
        gz="data/presabs/{feature}/{sample}.bed.gz",
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/summarize_coverage_over_feature/data/presabs/{feature}/{sample}.sum.tsv.benchmark.txt"
    log:
        "logs/summarize_coverage_over_feature/data/presabs/{feature}/{sample}.sum.tsv.log",
    threads: 1
    shell:
        """conifer-presabs summarize {input.gz} {output.tsv}"""
