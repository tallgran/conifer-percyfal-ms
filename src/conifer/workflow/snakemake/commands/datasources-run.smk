"""Data sources setup and provenance tracking.

Setup input data files from remote storage to local storage. There is
one generic rule that sets up input files as links or copies. The
setup behaviour is triggered by the file URI scheme, where file:/
triggers link creation and rsync:/ triggers file copy.

The datasources dictionary defines target:source mappings and is based
on input file names from specialized files (bamfiles, vcffiles) or
information in a datasources.yaml file.

"""
import re
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from urllib.parse import urlparse


##############################
# Configuration
##############################
configfile: Path("config/config.yaml")


envvars:
    "BAM_REMOTE_STORAGE",
    "VCF_REMOTE_STORAGE",
    "PACBIO_HIFI_STORAGE",
    "REFERENCE_STORAGE",
    "CHUNKED_REFERENCE_STORAGE",


BAM_REMOTE_STORAGE = os.environ["BAM_REMOTE_STORAGE"]
VCF_REMOTE_STORAGE = os.environ["VCF_REMOTE_STORAGE"]
datasources_envvars = {
    "$PACBIO_HIFI_STORAGE": os.environ["PACBIO_HIFI_STORAGE"],
    "$REFERENCE_STORAGE": os.environ["REFERENCE_STORAGE"],
    "$CHUNKED_REFERENCE_STORAGE": os.environ["CHUNKED_REFERENCE_STORAGE"],
}


# BAM files
BAM_OUT_PREFIX = "data/bam"
bamsources = {}
bamfiles = pd.read_table("resources/bamfiles.txt", header=None).iloc[:, 0].tolist()
for fn in bamfiles:
    bam = fn
    bai = fn + ".bai"
    bamsources[f"{BAM_OUT_PREFIX}/{bam}"] = f"{BAM_REMOTE_STORAGE}/{bam}"
    bamsources[f"{BAM_OUT_PREFIX}/{bai}"] = f"{BAM_REMOTE_STORAGE}/{bai}"


# VCF files
VCF_OUT_PREFIX = "data/vcf"
vcfsources = {}
vcffiles = pd.read_table("resources/vcffiles.txt", header=None).iloc[:, 0].tolist()
for fn in vcffiles:
    vcf = fn
    csi = fn + ".csi"
    vcfsources[f"{VCF_OUT_PREFIX}/{vcf}"] = f"{VCF_REMOTE_STORAGE}/{vcf}"
    vcfsources[f"{VCF_OUT_PREFIX}/{csi}"] = f"{VCF_REMOTE_STORAGE}/{csi}"

# datasources.yaml files
other = {}
if "datasources" in config:
    index = ["data"]
    with open(config["datasources"]) as fh:
        data = yaml.load(fh, yaml.Loader)
    for item in data:
        for key, value in datasources_envvars.items():
            if key in item["source"]:
                item["source"] = item["source"].replace(key, value)
    assert isinstance(data, list)
    df = pd.DataFrame(data)
    df.set_index(index, drop=False, inplace=True)
    df = df.replace({np.nan: None})
    df.index.names = index
    for target, source in zip(df.data, df.source):
        other[target] = [source]


# Make target:source mapping for all combinations
datasources = {**bamsources, **vcfsources, **other}


ALL = {
    "bamsources": bamsources.keys(),
    "vcfsources": vcfsources.keys(),
    "other": other.keys(),
}


rule all:
    input:
        **ALL,


rule all_bam:
    input:
        ALL["bamsources"],


rule all_vcf:
    input:
        ALL["vcfsources"],


rule all_other:
    input:
        ALL["other"],


##############################
# Atomic rules
##############################
rule setup_datasource:
    """Setup datasource from dictionary mapping. Generic rule for all datasources in datasources.yaml."""
    output:
        "{filename}",
    input:
        lambda wildcards: storage(datasources[wildcards.filename]),
    wildcard_constraints:
        filename="|".join(datasources.keys()),
    conda:
        "../envs/storage.yaml"
    benchmark:
        "benchmarks/setup_datasource/{filename}.benchmark.txt"
    log:
        "logs/setup_datasource/{filename}.log",
    threads: 1
    shell:
        """cp -d {input} {output}"""
