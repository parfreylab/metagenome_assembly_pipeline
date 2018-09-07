# Section 9: Functional annotation and basic plots

This next section uses two databases: the NCBI's Clusters of Orthologous Genes (COGs) and the Kyoto Encyclopedia of Genes and Genomes (KEGG) to functionally annotate contigs in each bin. After that, basic bar plots are generated for visualization of the annotations. Note however, that the `COGs` database is not as well maintained, with the latest update dating back to 2014.

## COGs database

Anvi'o provides a simple way to functionally annotate contigs by using NCBI's COGs. The first thing you need to do is to set up the COGs database. This can be done by running the `anvi-setup-ncbi-cogs` command. Note however, if you do not have superuser privileges, or you are working on the cluster, then you will need to pass in the `--cog-data-dir` flag:

```bash
anvi-setup-ncbi-cogs --cog-data-dir path/to/cogs/db
```

Once this is finished, annotate your contigs with `anvi-run-ncbi-cogs`:

```bash
anvi-run-ncbi-cogs -c my_contigs_db.db --cog-data-dir path/to/cogs/db
```

If you installed the COGs database in the default directory successfully, you do not need to pass in the `--cog-data-dir` flag. See the full guide [here][anvi-cog-annotation-link] for more detail.

## KEGG database

Since the KEGG database operates on a subscription based model, we need to work around this by using their online database (free), as well as their search tool, ghostKOALA (also free). There is a fairly comprehensive blog post [here][anvi-ghostkoala-kegg-annotation-link] on obtaining and importing KEGG annotations into Anvi'o. However, for this to work note the following:

- Rather than use `anvi-get-aa-sequences-for-gene-calls`, which no longer works, use:

```bash
anvi-get-sequences-for-gene-calls --get-aa-sequences -c contigs.db -o protein-sequences.fa
```

- Under the section ["generate the KEGG orthology table"][anvi-gen-kegg-ortho-table-link], the current `wget` command doesn't work (URL is a blank result). Instead, go directly to [this][kegg-ortho-table-link] link, and then on the top menu, click the option `Download htext`. You can continue following the guide with this file.

- Under the section ["parsing the results from GhostKOALA"][anvi-parse-ghostkoala-results-link], there is a line of code which is ran to add headers to the output file:

```bash
echo -e "contig\taccession_id" > .temp && cat user_ko.txt >> .temp && mv .temp user_ko.txt
```

This isn't necessary, and will actually cause an error when trying to import the annotations back to Anvi'o using `anvi-import-functions`. Specifically, refer to this [issue][ghost-koala-parser-issue-link] to see the error. As a result, you can skip that step.

## Plotting the gene calls

Once both COGs and KEGG functional annotations have been imported into Anvi'o, they can be exported bin by bin to make barplots. Custom scripts are available to make cursory barplots to visualize the functional annotations from each bin. These scripts are found under the [scripts][scripts-link] directory of this repository. Ultimately, this section goes from the gene calls files exported from Anvi'o to a tab-delimited text file required as input to the plotting script `barplot_rel_abundance.R`. To do this, gene calls from each bin need to be exported, and then the scripts can be run downstream for formatting and plotting.

If you need help with running the custom scripts, pass in the `-h` flag.

First, run `anvi-summarize` to get all gene calls from each bin.

```bash
anvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY-COG-KEGG -C BIN-COLLECTION-NAME
```

In the above example command, under the folder `SAMPLES-SUMMARY-COG-KEGG`, there will be another folder, `bin_by_bin`. This folder stores summary statistics for each individual bin, including the gene calls. Pull out the gene calls from bins of interest (or all of them) to a new folder to plot. If the names of the folders and files are left as default, then create a file `bin_names.txt`, with each line being the name of a bin. Example:

```bash
Bin_12
Bin_14
Bin_19
Bin_3
Bin_5
Bin_6
Bin_7
Refined_14
Refined_16
Refined_17
Refined_21
```

Run the loop to get all the gene calls files to a new directory:

```bash
# execute in folder SAMPLES-SUMMARY-COG-KEGG/bin_by_bin
mkdir -p gene_calls_files
while read name; do
    cp $name/$name-gene_calls.txt gene_calls_files/$name-gene_calls.txt
done < bin_names.txt
```

Take a look at one of the files. The header has 12 columns:

1. gene_callers_id	
2. contig	
3. start	
4. stop	
5. direction	
6. KeggGhostKoala	
7. KeggGhostKoala (ACCESSION)	
8. COG_FUNCTION	
9. COG_FUNCTION (ACCESSION)	
10. COG_CATEGORY	
11. COG_CATEGORY (ACCESSION)	
12. dna_sequence

The next plots will only use the `KeggGhostKoala (ACCESSION)` and `COG_CATEGORY (ACCESSION)` columns. Start with looking at the COG accessions. Each one is assigned a single letter, and this letter is in a category. To make sense of all these single letter assignments, they can be grouped together into higher categories, and these categories can be plotted for each bin to see the composition of the functional annotations.

Run the custom script `cog_anvio_parser.py` to parse the gene calls files from each bin. This will take in a directory of gene calls files, and output two tab-delimited .txt files with the COG categories on the left column, and every column following represents a bin. In each bin column, the numbers represent the number of functional annotations which are assigned to that category.

The first output file, `cogs_top_level_abundance.txt`, is grouping all of the annotations together under the broadest (highest level) COG categories. The second file, `cogs_low_level_abundance.txt` is for the lower level categories.

```bash
python cog_anvio_parser.py -d gene_calls_files/ -o out/
```

Next, notice that each entry under the `KeggGhostKoala (ACCESSION)` column is a KO identifier. Similar to the COGs parsing script, we can group these KO identifiers into higher categories as well using the script `group_kos.py`. This will produce a similar output to the previous script, but as KEGG has 4 hierarchical levels, 4 output tables can be generated. **NOTE however, that the default is to only output the first 3 levels.** This is because writing the fourth level has annotations at the protein level, which will result in too many colors in each bar to visualize. However, if you wish to make this table anyway, pass in the `-4` or `--include-fourth-level` flag.

The script `group_kos.py` takes in the KEGG Orthology csv file generated from the `KEGG-to-anvio` script, the directory containing all the gene calls files, and an optional output path.

```bash
python group_kos.py -k path/to/KeggOrthology_Table1.txt -g gene_calls_files/ -o out/
```

Since each bin will have an uneven amount of total counts per category, convert these values to percentages by using the script `convert_summary_table.py` on all of the output tables. This script can be run in a loop:

```bash
for f in *.txt; do
    python convert_summary_table.py -t $f -o out/
done
```

All the tables are now ready to be plotted. Run the Rscript `barplot_rel_abundance.R` with the converted tables as input. Note that the colors chosen for each category are **random**, so this may result in two categories sharing the same color. If this is the case, simply run the script multiple times until you are satisfied with the colors. Note that this script is also for a quick visualization of the annotations, and will need to be adapted for individual needs.

`barplot_rel_abundance.R` will make a barplot with each bar representing a bin (labelled with bin name) on the x axis. On the y axis will be the percentage of proteins in each category, marked by (hopefully) different colors.

`barplot_rel_abundance.R` takes in the file to be plotted, and optionally, titles for the graph, x and y axes labels, and and output directory.

```bash
Rscript barplot_rel_abundance.R -f summary_table_converted.txt -t GRAPH_TITLE -x X_AXIS_LABELS -y Y_AXIS_LABEL -o out/
```

Rerun this script for each level of the KEGG hierarchy you are interested in. If you are only interested in certain categories, then the tab delimited .txt file can be edited by removing rows you do not want to include in the plot. 

## Next step

Proceed to [section 10][section10-link].

[anvi-cog-annotation-link]: http://merenlab.org/2016/10/25/cog-annotation/
[anvi-ghostkoala-kegg-annotation-link]: http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
[anvi-gen-kegg-ortho-table-link]: http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/#generate-the-kegg-orthology-table
[kegg-ortho-table-link]: https://www.genome.jp/kegg-bin/get_htext?ko00001
[anvi-parse-ghostkoala-results-link]: http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/#parsing-the-results-from-ghostkoala
[ghost-koala-parser-link]: https://github.com/edgraham/GhostKoalaParser/blob/master/samples/KO_Orthology_ko00001.txt
[ghost-koala-parser-issue-link]: https://github.com/edgraham/GhostKoalaParser/issues/3
[scripts-link]: ../scripts
[section10-link]: ../section_10
