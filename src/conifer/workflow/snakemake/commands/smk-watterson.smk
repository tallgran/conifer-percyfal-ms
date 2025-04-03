import re
import sys
from conifer.snakemake.config import SchemaFiles
from conifer.snakemake import inputs
from conifer.snakemake.lmod import get_envmodules
from pathlib import Path
from snakemake.utils import validate
import pandas as pd
import numpy as np


envvars:
    "CONIFER_HOME",
    "PDC_TMP",
    "SNPEFF",
    "PABIES_FAI",


CONIFER_HOME = Path(os.environ["CONIFER_HOME"])
PDC_TMP = Path(os.environ["PDC_TMP"])
SNPEFF = Path(os.environ["SNPEFF"])
SNPEFF_DATA = Path(os.environ["SNPEFF_DATA"])
PABIES_FAI = Path(os.environ["PABIES_FAI"])

# Chunk size is 100 Mbp
OFFSET_BASES = 100_000_000


##############################
# Configuration
##############################
configfile: "config/config.yaml"


validate(config, schema=SchemaFiles.CONFIGURATION_SCHEMA)

samples = pd.read_table(config["samples"])
validate(samples, schema=SchemaFiles.SAMPLES_SCHEMA)


wildcard_constraints:
    data="data",
    interim="data/interim",
    reports="reports",
    sample=f"({'|'.join(samples.SM.tolist())})",


fai = pd.read_table(
    PABIES_FAI,
    header=None,
    names=["chrom", "length", "offset", "linebases", "linewidth"],
)
CHROMOSOMES = [x for x in fai["chrom"].tolist() if re.match(r"PA_chr\d+_\d+$", x)]

d4count = {
    "P.abies": "data/aggregate_coverage/P.abies/MQ10/count_ge3.p100.d4",
    "highcov": "data/aggregate_coverage/highcov/MQ10/count_ge3.p100.d4",
    "north": "data/aggregate_coverage/north/MQ10/count_ge3.p100.d4",
    "south": "data/aggregate_coverage/south/MQ10/count_ge3.p100.d4",
}
count_mask = {
    "P.abies": "data/mask/P.abies/count_ge3.filter_min528_max1070.p100.bed.gz",
    "highcov": "data/mask/highcov/count_ge3.filter_min148_max1070.p100.bed.gz",
    "north": "data/mask/north/count_ge3.filter_min13_max1070.p100.bed.gz",
    "south": "data/mask/south/count_ge3.filter_min13_max1070.p100.bed.gz",
}
sum_mask = {
    "P.abies": "data/mask/P.abies/sum.filter_min4248_max11258.p100.bed.gz",
    "highcov": "data/mask/highcov/sum.filter_min2000_max4500.p100.bed.gz",
    "north": "data/mask/north/sum.filter_min166_max512.p100.bed.gz",
    "south": "data/mask/south/sum.filter_min159_max514.p100.bed.gz",
}


# Define main ALL target
ALL = (
    [
        "data/watterson/P.abies/harmonic_count.bed.gz",
        "data/watterson/P.abies/mask.bed.gz",
    ]
    + [
        "data/watterson/highcov/harmonic_count.bed.gz",
        "data/watterson/highcov/mask.bed.gz",
    ]
    + [
        "data/watterson/P.abies/snpeff.genes.txt",
        "data/watterson/highcov/snpeff.genes.txt",
    ]
    + [
        "data/watterson/north/harmonic_count.bed.gz",
        "data/watterson/south/harmonic_count.bed.gz",
    ]
)


rule all:
    input:
        ALL,


##############################
# Atomic rules
##############################
rule merge_mask:
    """Merge count and sum mask files."""
    output:
        bed="data/watterson/{sampleset}/mask.bed.gz",
        gzi="data/watterson/{sampleset}/mask.bed.gz.gzi",
    input:
        count_mask=lambda wildcards: count_mask[wildcards.sampleset],
        sum_mask=lambda wildcards: sum_mask[wildcards.sampleset],
    envmodules:
        *get_envmodules(config, "bedtools"),
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/merge_mask/data/watterson/{sampleset}/mask.bed.gz.benchmark.txt"
    log:
        "logs/merge_mask/data/watterson/{sampleset}/mask.bed.gz.log",
    threads: 1
    shell:
        """
        bedtools intersect -a {input.count_mask} -b {input.sum_mask} | bgzip -c > {output.bed} | bgzip -i -o {output.bed}
        """


rule bcftools_subset:
    """Subset VCF to mask."""
    output:
        vcf=temp("data/watterson/{sampleset}/{chrom}.vcf.gz"),
        csi=temp("data/watterson/{sampleset}/{chrom}.vcf.gz.csi"),
    input:
        bed="data/watterson/{sampleset}/mask.bed.gz",
        vcf="data/vcf/{chrom}.biallelic.vcf.gz",
        fai="data/resources/Picab02_chromosomes.fasta.gz.fai",
        sampleset="resources/samplesets/samples-{sampleset}.keep.txt",
    params:
        offset=lambda wildcards: (
            int(re.search(r"_(\d+)$", wildcards.chrom).group(1)) - 1
        )
        * OFFSET_BASES,
        samples=lambda wildcards, input: ",".join(
            [
                x.capitalize()
                for x in pd.read_table(input.sampleset).SM.tolist()
                if not x.startswith("haploid")
            ]
        ),
    envmodules:
        *get_envmodules(config, "bcftools"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_subset/data/watterson/{sampleset}/{chrom}.vcf.gz.benchmark.txt"
    log:
        "logs/bcftools_subset/data/watterson/{sampleset}/{chrom}.vcf.gz.log",
    threads: 1
    shell:
        """bcftools view -i 'QUAL>=20' -s {params.samples} -c 1:minor -R {input.bed} {input.vcf} | \
        bcftools annotate -x 'FORMAT/PL' | \
        awk -v offset={params.offset} '{{if ($1 ~ /^#/) print; else {{gsub(/_[0-9]+$/, "", $1); printf("%s\\t%s\\t", $1, $2+offset); for (i=3; i<NF; i++) {{printf("%s\\t", $i);}} print $NF; }} }}' | \
        bcftools reheader -f {input.fai} | \
        bcftools view --write-index=csi -o {output.vcf} -O z 2>&1 > {log}
        """


rule snpeff:
    """Run snpeff on VCF."""
    output:
        vcf=temp("data/watterson/{sampleset}/{chrom}.snpeff.vcf.gz"),
        csi=temp("data/watterson/{sampleset}/{chrom}.snpeff.vcf.gz.csi"),
        csv="data/watterson/{sampleset}/{chrom}.snpeff.csv",
        txt="data/watterson/{sampleset}/{chrom}.snpeff.genes.txt",
    input:
        vcf="data/watterson/{sampleset}/{chrom}.vcf.gz",
        csi="data/watterson/{sampleset}/{chrom}.vcf.gz.csi",
    params:
        snpeff=SNPEFF,
    envmodules:
        *get_envmodules(config, "bcftools"),
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/snpeff/data/watterson/{sampleset}/{chrom}.snpeff.vcf.gz.benchmark.txt"
    log:
        "logs/snpeff/data/watterson/{sampleset}/{chrom}.snpeff.vcf.gz.log",
    threads: 16
    shell:
        """{params.snpeff} -Xmx16g -no-intergenic -no-downstream -no-upstream -csvStats {output.csv} pabi02 {input.vcf} | bcftools view --write-index=csi -O z -o {output.vcf} 2>&1 > {log}"""


rule harmonic_depth:
    """Calculate the mean harmonic depth."""
    output:
        bed="data/watterson/{sampleset}/harmonic_count.bed.gz",
        gzi="data/watterson/{sampleset}/harmonic_count.bed.gz.gzi",
    input:
        d4=lambda wildcards: d4count[wildcards.sampleset],
        cds="data/resources/regions/CDS.bed",
        count_mask=lambda wildcards: count_mask[wildcards.sampleset],
    params:
        mincov=lambda wildcards, input: re.search(
            r"min(\d+)_max(\d+)", Path(input.count_mask).stem
        ).group(1),
    envmodules:
        *get_envmodules(config, "d4tools"),
    conda:
        "../envs/d4tools.yaml"
    benchmark:
        "benchmarks/harmonic_depth/data/watterson/{sampleset}/harmonic_count.txt.benchmark.txt"
    log:
        "logs/harmonic_depth/data/watterson/{sampleset}/harmonic_count.txt.log",
    threads: 1
    shell:
        """d4tools view {input.d4} -R {input.cds} | awk '{{if ($4 >= {params.mincov}) print $0}}' | bgzip -i -o {output.bed} 2>&1 > {log}"""


rule concat_snpeff:
    """Concat snpeff resutls."""
    output:
        txt="data/watterson/{sampleset}/snpeff.genes.txt",
    input:
        txt=expand(
            "data/watterson/{{sampleset}}/{chrom}.snpeff.genes.txt", chrom=CHROMOSOMES
        ),
    benchmark:
        "benchmarks/concat_snpeff/data/watterson/{sampleset}/snpeff.genes.txt.benchmark.txt"
    log:
        "logs/concat_snpeff/data/watterson/{sampleset}/snpeff.genes.txt.log",
    threads: 1
    run:
        import pandas as pd

        dflist = [pd.read_table(f, sep="\t", skiprows=1) for f in input.txt]
        maxcol = max([len(df.columns) for df in dflist])
        i = min(np.where([len(df.columns) == maxcol for df in dflist])[0])
        data = pd.DataFrame(columns=dflist[i].columns)
        for i in range(len(dflist)):
            data = data.merge(dflist[i], how="outer")
        data.to_csv(output.txt, sep="\t", index=False)


rule fasta2axt:
    """Convert Fasta to AXT."""
    output:
        axt=temp("data/watterson/cds.axt"),
    input:
        SNPEFF_DATA / "cds.fa.gz",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/ fasta2axt/data/watterson/cds.axt.benchmark.txt"
    log:
        "logs/ fasta2axt/data/watterson/cds.axt.log",
    threads: 1
    shell:
        """conifer-fasta2axt {input} > {output.axt}"""


rule kaks:
    """Convert axt to kaks."""
    output:
        kaks="data/watterson/cds.axt.kaks",
    input:
        axt="data/watterson/cds.axt",
    envmodules:
        *get_envmodules(config, "conifer"),
    conda:
        "../envs/conifer.yaml"
    benchmark:
        "benchmarks/kaks/data/watterson/cds.axt.kaks.benchmark.txt"
    log:
        "logs/kaks/data/watterson/cds.axt.kaks.log",
    threads: 1
    shell:
        """KaKs -i {input.axt} -o {output.kaks} -m YN"""


localrules:
    concat_snpeff,
