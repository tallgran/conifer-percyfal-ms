rule bedtools_chromosome_genome_windows:
    """Make genome windows on chromosomes with bedtools"""
    output:
        bed="data/resources/windows/{prefix}.fa.fai.{window_size}.chromosomes.bed",
        names="data/resources/windows/{prefix}.fa.fai.{window_size}.chromosomes.bed.names",
    input:
        bed="data/resources/{prefix}.fa.fai.chromosomes.bed",
    conda:
        "../envs/bedtools.yaml"
    envmodules:
        *config["envmodules"]["bedtools"],
    benchmark:
        "benchmarks/bedtools_genome_windows/data/resources/windows/{prefix}.fa.fai.{window_size}.chromosomes.bed.benchmark.txt"
    log:
        "logs/bedtools_genome_windows/data/resources/windows/{prefix}.fa.fai.{window_size}.chromosomes.bed.log",
    threads: 1
    shell:
        """bedtools makewindows -w {wildcards.window_size} -b {input.bed} > {output.bed};
        cut -f 1 {output.bed} | uniq > {output.names}
        """


localrules:
    bedtools_chromosome_genome_windows,
