from conifer.snakemake.lmod import get_envmodules


rule bgzip_unzip_bed:
    """Unzip a bgziped bed file."""
    output:
        "{prefix}.bed",
    input:
        "{prefix}.bed.gz",
    envmodules:
        *get_envmodules(config, "tabix_bgzip"),
    conda:
        "../envs/tabixbgzip.yaml"
    benchmark:
        "benchmarks/bgzip_unzip/{prefix}.bed.benchmark.txt"
    log:
        "logs/bgzip_unzip/{prefix}.bed.log",
    threads: 1
    shell:
        """bgzip -d -c {input} > {output}"""


rule bgzip_bed:
    """bgzip and index a bed file."""
    output:
        bed="{prefix}.bed.gz",
        tbi="{prefix}.bed.gz.tbi",
    input:
        "{prefix}.bed",
    envmodules:
        *get_envmodules(config, "tabix_bgzip"),
    conda:
        "../envs/tabixbgzip.yaml"
    benchmark:
        "benchmarks/bgzip/{prefix}.bed.benchmark.txt"
    log:
        "logs/bgzip/{prefix}.bed.log",
    threads: 1
    shell:
        """bgzip -c {input} > {output.bed}
        tabix -p bed {output.bed}"""


localrules:
    bgzip_unzip_bed,
    bgzip_bed,
