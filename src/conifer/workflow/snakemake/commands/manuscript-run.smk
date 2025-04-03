"""Manuscript workflow.

Generate source data for the manuscript.
"""
import os
import re
import sys
from conifer.snakemake.config import SchemaFiles
from conifer.snakemake.lmod import get_envmodules
from pathlib import Path
from snakemake.utils import validate
import pandas as pd


envvars:
    "REFERENCE_DIR",
    "MANUSCRIPT_ROOT",


REFERENCE_DIR = Path(os.environ["REFERENCE_DIR"])
MANUSCRIPT_ROOT = Path(os.environ["MANUSCRIPT_ROOT"])
MAINGENOME_2024 = MANUSCRIPT_ROOT / "mainGenome2024"

##############################
# Input
##############################
PICAB = REFERENCE_DIR / "Picea-abies/v2.0"
PINSY = REFERENCE_DIR / "Pinus-sylvestris/v1.0"

picab = "picab02"
pinsy = "pinsy01"
species = [picab, pinsy]
window_size = [50000, 100000, 1000000, 10000000]
samplesets = [
    "P.abies",
    "highcov",
    "north",
    "south",
    # "recombination",
]
features = ["CDS", "UTR", "intron", "TE.gene", "pseudogene", "intergenic", "genome"]

INPUT = {
    "fai": {
        picab: PICAB / "fasta/Picab02_chromosomes.fasta.gz.fai",
        pinsy: PINSY / "fasta/Picab02_chromosomes.fasta.gz.fai",
    }
}


##############################
# Configuration
##############################
configfile: Path("config/config.yaml")


validate(config, schema=SchemaFiles.CONFIGURATION_SCHEMA)


wildcard_constraints:
    data="data",
    feature="|".join(features),
    interim="data/interim",
    reports="reports",
    sampleset="|".join(samplesets),
    species="|".join(species),
    window_size="|".join(map(str, window_size)),


# Define main ALL target
ALL = {
    "tbl-samplesets": MAINGENOME_2024 / "source_data/tbl-samplesets.csv",
    "tbl-feature-size": MAINGENOME_2024 / "source_data/tbl-feature-size.csv",
    "tbl-diversity-summary": MAINGENOME_2024 / "source_data/tbl-diversity-summary.csv",
    "tbl-pabies-diversity-summary-twocol": MAINGENOME_2024
    / "source_data/tbl-pabies-diversity-summary-twocol.csv",
    "tbl-highcov-diversity-summary-twocol": MAINGENOME_2024
    / "source_data/tbl-highcov-diversity-summary-twocol.csv",
    "tbl-watterson-snpeff": MAINGENOME_2024 / "source_data/tbl-watterson-snpeff.csv",
    "fig-individual-coverage": expand(
        MAINGENOME_2024 / "source_data/fig-individual-coverage-{sampleset}.csv",
        sampleset=["highcov", "north", "south"],
    ),
    "nucleotide-diversity": expand(
        MAINGENOME_2024
        / "figshare/genome-scale-variation/nucleotide-diversity-{sampleset}-w{window_size}.tsv.gz",
        sampleset=["P.abies", "highcov", "north", "south"],
        window_size=[50000, 100000],
    ),
    "fig-nucleotide-diversity": expand(
        MAINGENOME_2024
        / "source_data/fig-nucleotide-diversity-{sampleset}-w{window_size}-0.5-1000.tsv.gz",
        sampleset=["P.abies", "highcov", "north", "south"],
        window_size=[50000, 100000],
    ),
    "fig-genic-proximity": expand(
        MAINGENOME_2024
        / "source_data/fig-diversity-genic-proximity-{sampleset}.csv.gz",
        sampleset=["highcov"],
    ),
}


rule all:
    input:
        **ALL,


##############################
# Atomic rules
##############################
rule make_windows:
    """Make windows from input fasta index file"""
    output:
        bed="results/manuscript/reference/{species}.{window_size}.bed",
    input:
        fai=lambda wildcards: INPUT["fai"][wildcards.species],
    params:
        cat=lambda wildcards, input: "zcat" if input.fai.endswith(".gz") else "cat",
    envmodules:
        *get_envmodules(config, "bedtools"),
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/make_windows/results/manuscript/reference/{species}.{window_size}.bed.benchmark.txt"
    log:
        "logs/make_windows/results/manuscript/reference/{species}.{window_size}.bed.log",
    threads: 1
    shell:
        """{params.cat} {input.fai} | awk '{{print $1"\\t"$2}}' | \
        bedtools makewindows -w {wildcards.window_size} -g /dev/stdin > \
        {output.bed} > {log} 2>&1
        """


rule samplesets_table:
    """Make samplesets table"""
    output:
        MAINGENOME_2024 / "source_data/tbl-samplesets.csv",
    input:
        expand("resources/samplesets/samples-{sampleset}.tsv", sampleset=samplesets),
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/samplesets_table/source_data/tbl-samplesets.csv.benchmark.txt"
    log:
        "logs/samplesets_table/source_data/tbl-samplesets.csv.log",
    threads: 1
    shell:
        """conifer-manuscript-source-data make-samplesets-table"""


rule make_samplesheet_table:
    """Compile all samples with read counts, alignment rates, and coverage"""
    output:
        csv=MAINGENOME_2024 / "source_data/tbl-samplesheet.csv",
    input:
        tsv="resources/samplesheet.tsv",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/make_samplesheet_table/source_data/tbl-samplesheet.csv.benchmark.txt"
    log:
        "logs/make_samplesheet_table/source_data/tbl-samplesheet.csv.log",
    threads: 1
    shell:
        """echo {output.csv}"""


rule feature_size:
    """Make feature size output"""
    output:
        MAINGENOME_2024 / "source_data/tbl-feature-size.csv",
        MAINGENOME_2024 / "source_data/tbl-diversity-summary.csv",
        MAINGENOME_2024 / "source_data/tbl-pabies-diversity-summary-twocol.csv",
        MAINGENOME_2024 / "source_data/tbl-highcov-diversity-summary-twocol.csv",
    input:
        expand(
            "results/diversity/MQ10/{sampleset}/sites.pi.w1000000.summary.csv.gz",
            sampleset=samplesets[0:4],
        ),
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/feature_size/source_data/tbl-feature-size.csv.benchmark.txt"
    log:
        "logs/feature_size/source_data/tbl-feature-size.csv.log",
    threads: 1
    shell:
        """conifer-manuscript-source-data compile-nucleotide-diversity"""


rule watterson_snpeff:
    """Make Watterson output based on snpEff results"""
    output:
        csv=MAINGENOME_2024 / "source_data/tbl-watterson-snpeff.csv",
    input:
        pabies="data/watterson/P.abies/snpeff.genes.txt",
        highcov="data/watterson/highcov/snpeff.genes.txt",
        kaks="data/watterson/cds.axt.kaks",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/watterson_snpeff/output.benchmark.txt"
    log:
        "logs/watterson_snpeff/output.log",
    threads: 1
    shell:
        """
        conifer-watterson snpeff -H  data/watterson/P.abies/snpeff.genes.txt data/watterson/cds.axt.kaks -d P.abies -n 1529 -p 0.726137745694394 > {output.csv};
        conifer-watterson snpeff data/watterson/highcov/snpeff.genes.txt data/watterson/cds.axt.kaks -d highcov -n 567 -p 0.7020538909933343 >> {output.csv}
        """


rule individual_coverage:
    """Compile individual coverage data"""
    output:
        MAINGENOME_2024 / "source_data/fig-individual-coverage-{sampleset}.csv",
    input:
        expand(
            "data/aggregate_coverage/{{sampleset}}/MQ10/count_ge3.p100.d4.{feature}.hist",
            feature=features,
        ),
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/individual_coverage/source_data/fig-individual-coverage-{sampleset}.csv.benchmark.txt"
    log:
        "logs/individual_coverage/source_data/fig-individual-coverage-{sampleset}.csv.log",
    threads: 1
    shell:
        """conifer-manuscript-source-data compile-coverage --samplesets {wildcards.sampleset}"""


rule nucleotide_diversity:
    """Genome scale nucleotide diversity"""
    output:
        MAINGENOME_2024
        / "figshare/genome-scale-variation/nucleotide-diversity-{sampleset}-w{window_size}.tsv.gz",
    input:
        "results/diversity/MQ10/{sampleset}/sites.pi.w{window_size}.summary.csv.gz",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/nucleotide_diversity/figshare/genome-scale-variation/nucleotide-diversity-{sampleset}-w{window_size}.tsv.gz.benchmark.txt"
    log:
        "logs/nucleotide_diversity/figshare/genome-scale-variation/nucleotide-diversity-{sampleset}-w{window_size}.tsv.gz.log",
    threads: 1
    shell:
        """conifer-manuscript-source-data compile-genome-scale-nucleotide-diversity --samplesets {wildcards.sampleset} --window-size {wildcards.window_size}"""


rule fig_nucleotide_diversity:
    """Generate figure nucleotide diversity data"""
    output:
        MAINGENOME_2024
        / "source_data/fig-nucleotide-diversity-{sampleset}-w{window_size}-0.5-1000.tsv.gz",
    input:
        "results/diversity/MQ10/{sampleset}/sites.pi.w{window_size}.summary.csv.gz",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/fig_nucleotide_diversity/source_data/fig-nucleotide-diversity-{sampleset}-w{window_size}-0.5-1000.tsv.gz.benchmark.txt"
    log:
        "logs/fig_nucleotide_diversity/source_data/fig-nucleotide-diversity-{sampleset}-w{window_size}-0.5-1000.tsv.gz.log",
    threads: 1
    shell:
        """conifer-manuscript-source-data compile-filtered-nucleotide-diversity --samplesets {wildcards.sampleset} --window-size {wildcards.window_size}"""


rule fig_diversity_genic_proximity:
    """Generate figure diversity genic proximity data"""
    output:
        MAINGENOME_2024 / "source_data/fig-diversity-genic-proximity-{sampleset}.csv.gz",
    input:
        cds="results/diversity/MQ10/{sampleset}/sites.pi.w1000.e200000.CDS.summary.csv.gz",
        gene="results/diversity/MQ10/{sampleset}/sites.pi.w1000.e200000.gene.summary.csv.gz",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/fig_nucleotide_diversity/source_data/fig-diversity-genic-proximity-{sampleset}.csv.gz.benchmark.txt"
    log:
        "logs/fig_nucleotide_diversity/source_data/fig-diversity-genic-diversity-{sampleset}.csv.gz.log",
    threads: 1
    shell:
        """conifer-manuscript-source-data compile-genic-proximity --samplesets {wildcards.sampleset}"""
