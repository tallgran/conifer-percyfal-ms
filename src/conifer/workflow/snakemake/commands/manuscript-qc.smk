import re
import sys
from conifer.snakemake.config import SchemaFiles
from pathlib import Path
from snakemake.utils import validate
import pandas as pd


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


# Define main ALL target
ALL = {
    #    'qualimap': expand("results/manuscript/qc/qualimap/{sample}.html", sample=samples.SM.tolist()[-1]),
    #'fastqc': expand("results/manuscript/qc/fastqc/{sample}.html", sample=samples.SM.tolist()[-1]),
    "flagstat": expand(
        "results/manuscript/qc/samtools_flagstat/{sample}.txt",
        sample=samples.SM.tolist()[-1],
    ),
}


rule all:
    input:
        **ALL,


##############################
# Atomic rules
##############################
rule qualimap:
    """Collect qualimap results"""
    output:
        html="results/manuscript/qc/qualimap/{sample}.html",
        pdf="results/manuscript/qc/qualimap/{sample}.pdf",
        txt="results/manuscript/qc/qualimap/{sample}.txt",
        zip="results/manuscript/qc/qualimap/{sample}.zip",
    input:
        bam="data/bam/dedup.ir.{sample}.bam",
        bai="data/bam/dedup.ir.{sample}.bam.bai",
    benchmark:
        "benchmarks/qualimap/results/manuscript/qc/qualimap/{sample}.txt"
    log:
        "logs/qualimap/results/manuscript/qc/qualimap/{sample}.log",
    threads: 256
    shell:
        """qualimap bamqc -nt {threads} -bam {input.bam} -outdir results/manuscript/qc/qualimap/{wildcards.sample}"""


rule fastqc:
    """Run FastQC on raw data"""
    output:
        html="results/manuscript/qc/fastqc/{sample}.html",
    input:
        bam="data/bam/dedup.ir.{sample}.bam",
        bai="data/bam/dedup.ir.{sample}.bam.bai",
    benchmark:
        "benchmarks/fastqc/results/manuscript/qc/fastqc/{sample}.benchmark.txt"
    log:
        "logs/fastqc/results/manuscript/qc/fastqc/{sample}.log",
    threads: 256
    shell:
        """samtools view -h {input.bam} | fastqc -o results/manuscript/qc/fastqc/ -t {threads} stdin"""


rule samtools_flagstat:
    """Run samtools flagstat on bam file"""
    output:
        txt="results/manuscript/qc/samtools_flagstat/{sample}.txt",
    input:
        bam="data/bam/dedup.ir.{sample}.bam",
        bai="data/bam/dedup.ir.{sample}.bam.bai",
    benchmark:
        "benchmarks/samtools_flagstat/results/manuscript/qc/samtools_flagstat/{sample}.txt.benchmark.txt"
    log:
        "logs/samtools_flagstat/results/manuscript/qc/samtools_flagstat/{sample}.txt.log",
    threads: 256
    shell:
        """samtools flagstat -@ {threads} {input.bam} > {output.txt}"""
