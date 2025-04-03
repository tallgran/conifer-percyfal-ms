Coverage histogram
==================

{% if snakemake.wildcards.prefix.startswith("count") %}

Summary of number of samples with **at least** {{
snakemake.params.cutoff }}X coverage at each position (line plot, left
y-axis). The height of the y-axis displays the number of bases that
are covered by X samples (x-axis).

The data is also binned (red bar chart, right y-axis). Each bin sums
the number of bases that are covered by X samples.

{% elif snakemake.wildcards.prefix.startswith("sum") %}

Summary of the number of bases (first y-axis) with a given coverage
(x-axis). The plot has been truncated at the 99th upper percentile.
The histogram shows binned data, where each bin sums the number of
bases (second y-axis) within a given sequence coverage (x-axis).

This plot can be used to generate an accessibility mask based on a
coverage threshold, such that regions that lie outside the threshold
are either due to poor sample coverage or related to excessive mapping
to repetitive sequence.

{% else %}

{% endif %}
