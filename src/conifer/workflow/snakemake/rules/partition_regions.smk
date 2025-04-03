import numpy as np
from conifer.snakemake import inputs


rule fai_to_bed:
    output:
        "{prefix}.fai.bed",
    input:
        "{prefix}.fai",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/{prefix}.fai.bed.log",
    shell:
        "conifer-fai-to-bed {input} {output}"


rule partition_bed:
    output:
        "{prefix}.{p}-of-{partitions}.{method}.bed",
    input:
        "{prefix}.bed",
    wildcard_constraints:
        p="[0-9]+",
        partitions="[0-9]+",
        method="(sequential|greedy)",
    params:
        options=lambda wildcards: ""
        if (wildcards.method == "sequential")
        else "--greedy",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/{prefix}.{p}-of-{partitions}.{method}.bed",
    shell:
        "conifer-partition-bed {params.options} -p {wildcards.p} -n {wildcards.partitions} -o {output} {input}"


rule partition_samples_on_populationtype:
    output:
        "resources/samplesets/samples-{PopulationType}-{i}-of-{n}.tsv",
    input:
        "resources/samplesheet.tsv",
    wildcard_constraints:
        i="[0-9]+",
        n="[0-9]+",
        PopulationType="(RangeWide|BreedingPopulation|Outgroup)",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/resources/samplesets/samples-{PopulationType}-{i}-of-{n}.tsv.log",
    shell:
        "conifer-partition-samples -p 108 -w {wildcards.i} -o resources/samplesets -c PopulationType -i PopulationType:{wildcards.PopulationType} {input}"


rule partition_samples_on_populationtype_single:
    output:
        "resources/samplesets/samples-{PopulationType}.tsv",
    input:
        "resources/samplesheet.tsv",
    wildcard_constraints:
        PopulationType="(RangeWide|BreedingPopulation|Outgroup)",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/resources/samplesets/samples-{PopulationType}.tsv.log",
    shell:
        "conifer-partition-samples -s -o resources/samplesets -c PopulationType -i PopulationType:{wildcards.PopulationType} {input}"


rule partition_samples_P_abies:
    """Output all P.abies samples except the haploid"""
    output:
        "resources/samplesets/samples-P.abies.tsv",
    input:
        "resources/samplesheet.tsv",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/resources/samplesets/samples-P.abies.tsv.log",
    shell:
        "conifer-partition-samples -c Species -s -O {output} -i Species:P.abies -e SM:haploid_ERX242653 {input}"


rule partition_samples_on_populationtype_all:
    input:
        expand(
            "resources/samplesets/samples-BreedingPopulation-{i}-of-7.tsv",
            i=np.arange(7) + 1,
        ),
        expand(
            "resources/samplesets/samples-RangeWide-{i}-of-3.tsv", i=np.arange(3) + 1
        ),
        expand("resources/samplesets/samples-Outgroup-1-of-1.tsv"),
        expand(
            "resources/samplesets/samples-{PopulationType}.tsv",
            PopulationType=["RangeWide", "BreedingPopulation", "Outgroup"],
        ),


rule partition_samples_on_coverage:
    output:
        "{resources}/samplesets/samples-{sampleset}-coverage-min{covmin}-max{covmax}.tsv",
    input:
        coverage=inputs.partition_samples_on_coverage_input,
        samplesheet="{resources}/samplesets/samples-{sampleset}.tsv",
    wildcard_constraints:
        covmin=r"[0-9]+\.[0-9]+",
        covmax=r"[0-9]+\.[0-9]+",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/{resources}/samplesets/samples-{sampleset}-coverage-min{covmin}-max{covmax}.tsv.log",
    shell:
        "conifer-partition-samples -s -o {wildcards.resources}/samplesets --covmin {wildcards.covmin} --covmax {wildcards.covmax} {input.samplesheet} {input.coverage}"


rule partition_samples_on_coverage_omit_samples:
    output:
        "resources/samplesets/samples-{label}-coverage-min{covmin}-max{covmax}.tsv",
    input:
        coverage=inputs.partition_samples_on_coverage_input,
        samplesheet="resources/samplesheet.tsv",
        omit="resources/{label}.txt",
    wildcard_constraints:
        covmin=r"[0-9]+\.[0-9]+",
        covmax=r"[0-9]+\.[0-9]+",
    conda:
        "../envs/conifer.yaml"
    log:
        "logs/resources/samplesets/samples-{label}-coverage-min{covmin}-max{covmax}.tsv.log",
    shell:
        "conifer-partition-samples -s --omit-samples {input.omit} -o resources/samplesets --covmin {wildcards.covmin} --covmax {wildcards.covmax} {input.samplesheet} {input.coverage} > {log} 2>&1"


rule copy_reference_bed:
    output:
        "data/resources/{npartitions}/pabies-2.0.fa.fai.bed",
    input:
        "data/resources/pabies-2.0.fa.fai.bed",
    log:
        "logs/data/resources/{npartitions}/pabies-2.0.fa.fai.bed.log",
    conda:
        "../envs/conifer.yaml"
    shell:
        "cp {input} {output}"


localrules:
    fai_to_bed,
    partition_bed,
    partition_samples_on_populationtype,
    partition_samples_on_populationtype_single,
    partition_samples_on_coverage,
    partition_samples_on_coverage_omit_samples,
    partition_samples_P_abies,
