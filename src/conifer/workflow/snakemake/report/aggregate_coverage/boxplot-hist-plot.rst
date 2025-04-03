{% if snakemake.wildcards.plot_type.startswith("boxplot") %}
Coverage boxplot
==================
{% elif snakemake.wildcards.plot_type.startswith("violin") %}
Coverage violin plot
====================
{% endif %}

{% if snakemake.wildcards.plot_type.startswith("boxplot") %}

Boxplots of coverage distributions over features where coverage is the
number of samples with **at least** {{ snakemake.wildcards.coverage }}X
coverage at each position.

{% elif snakemake.wildcards.plot_type.startswith("violin") %}

Violin plots of coverage distributions over features where coverage is the
number of samples with **at least** {{ snakemake.wildcards.coverage }}X
coverage at each position.

{% else %}

{% endif %}

The sample size refers to the number of subsampled data points that
constitute the categories.
