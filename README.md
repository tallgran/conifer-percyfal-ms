<!--
[![PyPI](https://img.shields.io/pypi/v/conifer.svg)](https://pypi.python.org/pypi/conifer)
-->
<!--
[![CI](https://github.com/NBISweden/conifer/actions/workflows/ci.yml/badge.svg)](https://github.com/NBISweden/conifer/actions/workflows/ci.yml)
-->
<!--
[![BioConda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/conifer/README.html)
-->
[![DOI](https://zenodo.org/badge/960077764.svg)](https://doi.org/10.5281/zenodo.15141800)


# conifer

Conifer population genomic analyses.

# Installation

This repo was initialized with the cookiecutter
[cookiecutter-nbis-project](https://github.com/percyfal/cookiecutter-nbis-project).
It is recommended to install in a local conda environment. The repo
has to be under source control, whereafter it can be installed:

    git init
    git add -f .
    python -m pip install -e .

This will install the python module `conifer`. Type

    conifer

for help on the analyses and how to run tools.

## Conda environment

Binary requirements are listed in `environment.yml` and can be
installed by creating a dedicated conda environment:

    conda env create -n conifer -f environment.yml

See next section for information specific to the UPPMAX HPC compute on
which analyses were run.

# Project setup

Before doing any analyses source data need to be setup. This is
accomplished by running some of the provided scripts.

## Sample sheet

Run

    conifer tools samplesheet

or

    conifer-samplesheet

to generate samplesheet `resources/samplesheet.tsv`.

## Datasources

Data provenance is tracked in `resources/datasources.yaml` and setup
using the [datasources](https://github.com/percyfal/datasources-smk)
workflow.

To setup datasources run

    conifer datasources

# Running analyses

Type `conifer --help` for more detailed documentation.

## 0a. Initialize sample sets

`conifer smk init-samplesets` will initialize and define sample sets
used for comparisons. See output in `resources`.

## 0b. Generate region and feature files

`conifer smk annotation` generates region and feature files for
downstream analyses. See output in `data/resources`.

## 1. Coverage

`conifer smk coverage` runs mosdepth and generates individual coverage
profiles and reports for all samples. See `data/mosdepth_coverage` for
output.

Note that it may take some time to generate a list of jobs since this
is run on all ~1000 samples.

`conifer smk coverage --report` will generate summary reports and sync
them to webexport.

## 2. Aggregate coverage

`conifer smk aggregate-coverage` aggregates and summarizes individual
coverage profiles. The output is used to (manually) determine coverage
thresholds for accessibility masks. The thresholds are entered in the
configuration file; if a sampleset has missing thresholds, the
diversity calculation will not run.

`conifer smk aggregate-coverage --report` will generate summary
reports and sync them to webexport.

## 3. Diversity stats

`conifer smk diversity` generates masks files and vcftools results and
combines these to generate normalized output in `results/diversity`.

`conifer smk diversity --report` will generate summary reports and sync
them to webexport.

## 4. Synonymous / non-synonymous Watterson

`conifer smk watterson` generates an alternative calculation of
Watterson's theta by running SnpEff on coding sequences to tabulate
sense / missense variants, followed by normalization of synonymous /
non-synonymous "space" as estimated with KaKs

## 5. Presence absence data

`conifer smk presabs` generates mosdpeth coverage histogram input data
for `notebooks/presabs.ipynb`.

# Contact

Questions and discussion regarding these analyses may be directed to
<per.unneberg@scilifelab.se>
