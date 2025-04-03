"""Coverage module

## About

Collect coverage and mapping information for each sample. The workflow
runs mosdepth and samtools on all samples for different
MAPPING_QUALITY cutoffs.

The results are output in d4-format. The d4-format repository ships
the script d4tools which is used to view coverage files and generate
histogram data. The result files are used to generate plots for the
snakemake report.

Note that it may take some time to generate a list of all jobs as
there are

## Results

### mosdepth

Coverage results are output in
data/mosdepth_coverage/MQ{MAPPING_QUALITY}.

### samtools stats

Mapping statistics are output in results/samtools_stats/MQ{MAPPING_QUALITY}.

"""
import re
import sys
import math
import socket
from conifer.snakemake.config import SchemaFiles
from conifer.snakemake import inputs
from conifer.snakemake.lmod import get_envmodules
from pathlib import Path
from snakemake.utils import validate
import pandas as pd
import numpy as np
import pypandoc


##############################
# Configuration
##############################
configfile: Path("config/config.yaml")


config["__doc__"] = pypandoc.convert_text(__doc__, "rst", format="md")


report: "../report/coverage/workflow.rst"


validate(config, schema=SchemaFiles.CONFIGURATION_SCHEMA)

samples = pd.read_table(config["samples"])
validate(samples, schema=SchemaFiles.SAMPLES_SCHEMA)

REFERENCE = Path(config["ref"])


wildcard_constraints:
    data="data",
    interim="data/interim",
    mapping_quality="[0-9]+",
    partition_method="(sequential|greedy)",
    partition="[0-9]+",
    npartitions="[0-9]+",
    samplepartition="[0-9]+",
    min_quality="[0-9]+",
    reports="reports",
    results="results",
    sample=f"({'|'.join(samples.SM.tolist())})",
    sample_size="[0-9]+",


ALL = dict(
    config="results/config/coverage.config.yaml",
    multiqc=[],
    plots=[],
)
# Mapping analyses
for mq, conf in config["coverage"]["MQ"].items():
    ALL["multiqc"] += expand("reports/coverage/multiqc_MQ{MQ}.html", MQ=mq)
    if not "plots" in conf.keys():
        continue
    SM = conf["plots"].get("samples", samples.SM.tolist())
    ALL["plots"] += expand(
        "data/mosdepth_coverage/MQ{MQ}/{sample}.per-base.d4.hist.svg", MQ=mq, sample=SM
    )
    ALL["plots"] += expand(
        "data/mosdepth_coverage/MQ{MQ}/{sample}.per-base.d4.{window_size}.contigs.png",
        MQ=mq,
        window_size=conf["plots"]["window_size"],
        sample=SM,
    )
    SM = conf["plots"].get("html_samples", [])
    ALL["plots"] += expand(
        "data/mosdepth_coverage/MQ{MQ}/{sample}.per-base.d4.{window_size}.contigs.html",
        MQ=mq,
        window_size=conf["plots"]["window_size"],
        sample=SM,
    )
    ALL["plots"] += expand(
        "data/mosdepth_coverage/MQ{MQ}/{sample}.per-base.d4.hist.html", MQ=mq, sample=SM
    )


if config.get("__test__") is True:

    include: "test-coverage-config.smk"


rule all:
    input:
        **ALL,


##############################
# Atomic rules
##############################
rule mosdepth_d4:
    output:
        globaldist="{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.mosdepth.global.dist.txt",
        summary="{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.mosdepth.summary.txt",
        perbase="{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4",
    input:
        bam="{data}/bam/dedup.ir.{sample}.bam",
        bai="{data}/bam/dedup.ir.{sample}.bam.bai",
    params:
        prefix=lambda wildcards: f"{wildcards.data}/mosdepth_coverage/MQ{wildcards.mapping_quality}/{wildcards.sample}",
    threads: 10
    resources:
        mem_mb=60000,
    conda:
        "../envs/mosdepth.yaml"
    envmodules:
        *get_envmodules(config, "mosdepth"),
    log:
        "logs/{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.mosdepth.log",
    priority: 100
    benchmark:
        "logs/benchmark/{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.mosdepth.benchmark.txt"
    shell:
        "mosdepth --d4 -x {params.prefix} -Q {wildcards.mapping_quality} {input.bam} -t {threads}"


rule d4_stat:
    """Calculate statistics for d4 file"""
    output:
        "{data}/mosdepth_coverage/{prefix}.per-base{suffix}.d4.{stat}",
    input:
        "{data}/mosdepth_coverage/{prefix}.per-base{suffix}.d4",
    conda:
        "../envs/d4tools.yaml"
    wildcard_constraints:
        stat="(mean|median|hist)",
        suffix="(.sum|)",
    params:
        maxbin=lambda wildcards: 100000 if wildcards.suffix == ".sum" else 10000,
    log:
        "logs/{data}/mosdepth_coverage/{prefix}.per-base{suffix}.d4.{stat}.log",
    benchmark:
        "logs/benchmark/{data}/mosdepth_coverage/{prefix}.per-base{suffix}.d4.{stat}.benchmark.txt"
    threads: 4
    shell:
        "d4tools stat --max-bin {params.maxbin} -H -t {threads} -s {wildcards.stat} {input} > {output}"


rule multiqc_mosdepth_coverage_collect_inputs:
    """Collect mosdepth coverage input files for multiqc report"""
    output:
        txt="{data}/mosdepth_coverage/MQ{mapping_quality}/multiqc_input.txt",
    input:
        lambda wildcards: inputs.multiqc_coverage_input(
            wildcards,
            config["coverage"]["MQ"][int(wildcards.mapping_quality)].get(
                "samples", samples.SM
            ),
        ),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/multiqc_mosdepth_coverage/{data}/mosdepth_coverage/MQ{mapping_quality}/multiqc_input.txt.benchmark.txt"
    log:
        "logs/multiqc_mosdepth_coverage/{data}/mosdepth_coverage/MQ{mapping_quality}/multiqc_input.txt.log",
    threads: 1
    shell:
        """
        echo {input} | tr " " "\n" > {output.txt}
        """


rule multiqc:
    output:
        html="reports/coverage/multiqc_MQ{mapping_quality}.html",
    input:
        txt="data/mosdepth_coverage/MQ{mapping_quality}/multiqc_input.txt",
    envmodules:
        *get_envmodules(config, "multiqc"),
    params:
        configfile=inputs.multiqc_configfile,
    conda:
        "../envs/multiqc.yaml"
    benchmark:
        "benchmarks/multiqc/coverage/reports/multiqc_MQ{mapping_quality}.html.benchmark.txt"
    log:
        "logs/multiqc/coverage/reports/multiqc_MQ{mapping_quality}.html.log",
    threads: 1
    shell:
        """
        multiqc --interactive -f -l {input.txt} --outdir $(dirname {output.html}) -n $(basename {output.html}) > {log} 2>&1
        """


rule report_mosdepth_coverage_plot:
    """Make coverage plot for snakemake report"""
    output:
        report(
            "{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.hist.{ext}",
            caption="../report/coverage/coverage-plot.rst",
            category="MQ{mapping_quality}",
            subcategory="Coverage",
        ),
    input:
        "{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.hist",
    params:
        backend=lambda wildcards: {"svg": "mpl", "html": "plotly"}[wildcards.ext],
        max_bin="--max-bin 0.99",
        title=lambda wildcards: f"{wildcards.sample} coverage histogram",
        num_bins="--num-bins 20",
    conda:
        "../envs/plotting.yaml"
    benchmark:
        "benchmarks/mosdepth_coverage_plot/{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.hist.{ext}.benchmark.txt"
    log:
        "logs/mosdepth_coverage_plot/{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.hist.{ext}.log",
    threads: 1
    shell:
        """
        conifer-plot hist {input} -o {output} --ylab "Genome coverage (bp)" --title "{params.title}" --xlab "Depth of coverage (X)" --width 15 --height 6 --backend {params.backend} {params.max_bin} {params.num_bins}
        """


rule report_plot_contig_coverage:
    """Plot contig coverage in windows"""
    output:
        report(
            "{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.{window_size}.contigs.{ext}",
            caption="../report/coverage/coverage-contig-plot.rst",
            category="MQ{mapping_quality}",
            subcategory="Contig coverage",
        ),
    input:
        d4="{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4",
        bed=f"data/resources/windows/{REFERENCE.name}.fai.{{window_size}}.chromosomes.bed",
        names=f"data/resources/windows/{REFERENCE.name}.fai.{{window_size}}.chromosomes.bed.names",
    params:
        backend=lambda wildcards: {"png": "mpl", "html": "plotly"}[wildcards.ext],
    conda:
        "../envs/plotting.yaml"
    benchmark:
        "benchmarks/report_plot_contig_coverage/{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.{window_size}.contigs.{ext}.benchmark.txt"
    log:
        "logs/report_plot_contig_coverage/{data}/mosdepth_coverage/MQ{mapping_quality}/{sample}.per-base.d4.{window_size}.contigs.{ext}.log",
    threads: 1
    shell:
        """d4tools stat -t {threads} -r {input.bed} {input.d4} | bedtools sort -faidx {input.names} | conifer-plot coverage - {output} --title "{wildcards.sample}, window size {wildcards.window_size}" --width 15 --height 6 --backend {params.backend}"""


rule samtools_stats:
    """Collect samtools statistics"""
    output:
        "results/samtools_stats/{sample}.txt",
    input:
        bam="data/bam/dedup.ir.{sample}.bam",
        bai="data/bam/dedup.ir.{sample}.bam.bai",
        ref=REFERENCE,
    envmodules:
        *get_envmodules(config, "samtools"),
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools_stats/results/samtools_stats/{sample}.txt.benchmark.txt"
    log:
        "logs/samtools_stats/results/samtools_stats/{sample}.txt.log",
    threads: 1
    shell:
        """samtools stats -r {input.ref} {input.bam} > {output} 2> {log}"""


include: "../rules/partition_regions.smk"
include: "../rules/config.smk"
include: "../rules/windowing.smk"


localrules:
    copy_reference_bed,
    report_mosdepth_coverage_plot,


if config.get("__test__") is True:

    include: "test-coverage-setup.smk"
    include: "test-setup.smk"
