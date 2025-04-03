Coverage contig plot
=====================

{% set directory_string = snakemake.output | string() %}
{% set directory_name = directory_string.split('/') | last %}

{% if directory_name.startswith("count") %}

Coverage plot over all chromosomes where coverage is the number of
samples with **at least** {{ snakemake.params.coverage }}X coverage at each
position.

{% elif directory_name.startswith("sum") %}

Coverage plot over all chromosomes where coverage is the total
coverage for all samples at each position.

{% else %}

{% endif %}
