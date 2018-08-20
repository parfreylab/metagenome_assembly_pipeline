# Section 7: Binning quality evaluation

This section evaluates the quality of the binning approaches and bins obtained used in the previous section.

## Binning quality checking, lineage, and plots/statistics

We will use CheckM for quality check of bins. Use the [lineage workflow][checkm-lineage-wf-link] to place the genome bins into a reference genome tree and to identify marker genes and estimate contamination. Note that An'vio can also do completeness and redundancy estimates (via `anvi-summarize`), but these reports are likely to differ based on their different methods and different sets of marker genes.

Run the following command on the cluster, as CheckM has Linux dependencies. These example commands use bins from MetaBAT, but be sure to run these commands on the outputs of each of the binning programs you used in the previous section.

```bash
checkm lineage_wf -t 4 -x fa metabat_bins checkm/
```

Use CheckM for plots of bin quality. Specifically, use `bin_qa_plot` for a visual representation of completeness, contamination, and strain heterogeneity. [Link][checkm-bin-qa-plot-link] for description of plots.

```bash
checkm bin_qa_plot -x fa checkm metabat_bins plots
```

Use the CheckM utility command `coverage` to get coverage profiles for all sequences within the genome bins created. Coverage profiles are also required for a number of different plots produced by checkM.

```bash
checkm coverage -x fa metabat_bins coverage.tsv example_1.bam example_2.bam
```

Use the `profile` utility to produce a table with bin size, mapped reads, % mapped reads, % binned populations, and % community. The output defaults to `stdout`, so include the option `-f` (and `--tab_table`) to write to a file. The percentages indicate percentages of reads **mapped** to an assembly.

```bash
checkm profile coverage.tsv
```

## Next step

Proceed to [section 8][section8-link].

[checkm-lineage-wf-link]: https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow
[checkm-bin-qa-plot-link]: https://github.com/Ecogenomics/CheckM/wiki/Plots#bin_qa_plot
[section8-link]: ../section_8