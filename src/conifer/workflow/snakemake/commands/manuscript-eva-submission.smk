import os
import re
from pathlib import Path


envvars:
    "VCFDIR",
    "SCRATCH",
    "PICAB02_FA",
    "PICAB02_FAI",


# Chunk size is 100 Mbp
OFFSET_BASES = 100_000_000

VCFDIR = Path(os.environ["VCFDIR"])
SCRATCH = Path(os.environ["SCRATCH"])
WORKDIR = SCRATCH / "conifer" / "EVA"


workdir: WORKDIR


# Read the list of VCF files
vcf_files = list(VCFDIR.glob("*.vcf.gz"))
prefixes = [f.with_suffix("").stem for f in vcf_files]

print("===================")
print()
print("VCFDIR:", VCFDIR)
print("VCF file count:", len(vcf_files))
print("Output directory: ", os.path.abspath("."))
print()
print("===================")


wildcard_constraints:
    prefix="(" + "|".join(prefixes) + ")",


rule all:
    input:
        expand("{prefix}.vcf.gz.stats", prefix=prefixes),  # + ["all.vcf.validated.txt"]


##############################
# Atomic rules
##############################
rule copy_fasta_index:
    """Copy fasta index file to the working directory."""
    output:
        fai="Picab02_chromosomes.fasta.gz.fai",
    input:
        fai=os.environ["PICAB02_FAI"],
    benchmark:
        "benchmarks/link_fasta_index/Picab02_chromosomes.fasta.gz.fai.benchmark.txt"
    log:
        "logs/link_fasta_index/Picab02_chromosomes.fasta.gz.fai.log",
    threads: 1
    shell:
        """cp {input.fai} {output.fai}"""


rule make_unzipped_fasta_reference:
    """Zcat the fasta file and write to a new file."""
    output:
        fa="Picab02_chromosomes.fasta",
    input:
        fa=os.environ["PICAB02_FA"],
    benchmark:
        "benchmarks/make_unzipped_fasta_reference/Picab02_chromosomes.fasta.benchmark.txt"
    log:
        "logs/make_unzipped_fasta_reference/Picab02_chromosomes.fasta.log",
    threads: 1
    shell:
        """zcat {input.fa} > {output.fa}"""


rule convert_vcf_coordinates:
    """Add a multiple of OFFSET_BASES to the POS field and rename contig.
    """
    output:
        vcf="{prefix}.vcf.gz",
        csi="{prefix}.vcf.gz.csi",
    input:
        vcf=VCFDIR / "{prefix}.vcf.gz",
        fai="Picab02_chromosomes.fasta.gz.fai",
    envmodules:
        "bcftools/1.20",
    params:
        offset=lambda wildcards: (
            int(re.search(r"_(\d+)$", wildcards.prefix).group(1)) - 1
        )
        * OFFSET_BASES,
    benchmark:
        "benchmarks/convert_vcf_coordinates/{prefix}.vcf.gz.benchmark.txt"
    log:
        "logs/convert_vcf_coordinates/{prefix}.vcf.gz.log",
    threads: 1
    shell:
        """bcftools view {input.vcf} | grep -v "##bcftools\|reference=\|##ALT=<ID" | \
        bcftools annotate -x 'FORMAT/PL' --no-version | \
        awk -v offset={params.offset} '{{if ($1 ~ /^#/) print; else {{gsub(/_[0-9]+$/, "", $1); printf("%s\\t%s\\t", $1, $2+offset); for (i=3; i<NF; i++) {{printf("%s\\t", $i);}} print $NF; }} }}' | \
        bcftools reheader -f {input.fai} | \
        bcftools view --write-index=csi -o {output.vcf} -O z 2>&1 > {log}
        """


rule bcftools_stats:
    """Run bcftools stats on the VCF file."""
    output:
        stats="{prefix}.vcf.gz.stats",
    input:
        vcf="{prefix}.vcf.gz",
    envmodules:
        "bcftools/1.20",
    benchmark:
        "benchmarks/bcftools_stats/{prefix}.vcf.gz.stats.benchmark.txt"
    log:
        "logs/bcftools_stats/{prefix}.vcf.gz.stats.log",
    threads: 1
    shell:
        """bcftools stats {input.vcf} > {output.stats}
        """


# rule bcftools_annotate_and_rehead:
#     """Remove PL field from FORMAT and add INFO field from a file.
#     Rehead the VCF using the genome fasta index file. """
#     output:
#         vcf="annotate/{prefix}.vcf.gz",
#         csi="annotate/{prefix}.vcf.gz.csi",
#     input:
#         vcf="{prefix}.vcf.gz",
#         fai="Picab02_chromosomes.fasta.gz.fai",
#     benchmark: "benchmarks/bcftools_annotate_and_rehead/annotate/{prefix}.vcf.gz.benchmark.txt",
#     log: "logs/bcftools_annotate/annotate_and_rehead/{prefix}.vcf.gz.log",
#     threads: 1
#     shell:
#         """bcftools view {input.vcf} | grep -v "##bcftools\|reference=\|##ALT=<ID" | \
#         bcftools reheader -f {input.fai} | \
#         bcftools annotate -x 'FORMAT/PL' --no-version --write-index=csi -o {output.vcf} -O z 2>&1 > {log}
#         """


# rule vcf_validate:
#     """Validate VCF file with EVA vcf_validator_linux tool."""
#     output:
#         validated="{prefix}.vcf.gz.validated",
#     input:
#         vcf="{prefix}.vcf.gz",
#     benchmark: "benchmarks/vcf_validate/{prefix}.vcf.gz.validated.benchmark.txt",
#     log: "logs/vcf_validate/{prefix}.vcf.gz.validated.log",
#     threads: 1
#     shell:
#         """vcf_validator_linux -i {input.vcf} -o vcf_validation 2>&1 > {log};
#         touch {output.validated}
#         """


rule eva_sub_cli_validate:
    """Validate VCF file with EVA vcf_validator_linux tool.

    NB: has to be run manually because the tool requires a local conda
    environment.

    """
    output:
        txt="all.vcf.validated.txt",
    input:
        vcf=expand("{prefix}.vcf.gz", prefix=prefixes),
        xlsx="EVA_Submission_Spruce.xlsx",
        fa="Picab02_chromosomes.fasta",
        nfconf="nextflow.config",
    benchmark:
        "benchmarks/eva_sub_cli_validate/all.vcf.validated.benchmark.txt"
    log:
        "logs/eva_sub_cli_validate/all.vcf.validated.log",
    threads: 1
    shell:
        """eva-sub-cli.py --submission_dir submission --tasks validate --shallow --metadata_xlsx {input.xlsx} --vcf_files {input.vcf} --reference_fasta {input.fa} --nextflow-config {input.nfconf} 2>&1 > {log};"""
