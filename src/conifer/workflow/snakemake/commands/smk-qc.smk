import re
import sys
from conifer.snakemake.config import SchemaFiles
from conifer.snakemake.lmod import get_envmodules
from conifer.snakemake.inputs import multiqc_qc_input
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
ALL = [
    "results/qc/qualimap/pacbio-hifi/diploid_sorted_stats/genome_results.txt",
    "results/qc/bamtools/dedup.ir.diploid.stats.txt",
    "results/qc/bamtools/pacbio-hifi/diploid_sorted.stats.txt",
    "results/qc/multiqc/multiqc_report.html",
]  # + expand("results/qc/qualimap/{sample}_stats/genome_results.txt", sample=samples.SM.tolist())


rule all:
    input:
        ALL,


##############################
# Atomic rules
##############################
rule qualimap_bamqc:
    """Run qualimap bamqc on bam files"""
    output:
        txt="results/qc/qualimap/{prefix}_stats/genome_results.txt",
    input:
        bam="data/bam/{prefix}.bam",
    envmodules:
        *get_envmodules(config, "qualimap"),
    conda:
        "../envs/qualimap.yaml"
    benchmark:
        "benchmarks/qualimap_bamqc/results/qc/{prefix}_stats/genome_results.txt.benchmark.txt"
    resources:
        mem_mb=128000,
    log:
        "logs/qualimap_bamqc/results/qc/{prefix}_stats/genome_results.txt.log",
    threads: 20
    shell:
        """
        unset DISPLAY; qualimap bamqc --java-mem-size={resources.mem_mb}M -c -nt {threads} -bam {input.bam} -outdir results/qc/qualimap/{wildcards.prefix}_stats > {log} 2>&1
        """


rule bamtools_stats:
    """Run bamtools stats on bam file"""
    output:
        txt="results/qc/bamtools/{prefix}.stats.txt",
    input:
        bam="data/bam/{prefix}.bam",
    envmodules:
        *get_envmodules(config, "bamtools"),
    conda:
        "../envs/bamtools.yaml"
    benchmark:
        "benchmarks/bamtools_stats/results/qc/bamtools/{prefix}.stats.txt.benchmark.txt"
    log:
        "logs/bamtools_stats/results/qc/bamtools/{prefix}.stats.txt.log",
    threads: 1
    shell:
        """
        bamtools stats -in {input.bam} > {output.txt}
        """


rule multiqc:
    """Make multiqc report"""
    output:
        fof="results/qc/multiqc/fof.txt",
        html="results/qc/multiqc/multiqc_report.html",
    input:
        lambda wildcards: multiqc_qc_input(config),
    envmodules:
        *get_envmodules(config, "multiqc"),
    conda:
        "../envs/multiqc.yaml"
    benchmark:
        "benchmarks/multiqc/qc/results/multiqc_report.html.benchmark.txt"
    log:
        "logs/multiqc/qc/results/multiqc_report.html.log",
    threads: 1
    shell:
        """
        echo > {output.fof};
        for fn in {input}; do cat $fn >> {output.fof}; done;
        multiqc --interactive -f -l {output.fof} --outdir results/multiqc > {log} 2>&1
        """
