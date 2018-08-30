# Section 8: Binning, taxonomy and bin refinement in Anvi'o

This section imports contigs generated from previous sections into Anvi'o, bins with CONCOCT, assigns taxonomy to bins, refines bins, and compares multiple binning approaches. Since the majority of these steps already have good documentation on the Meren lab website, this section will mostly link to those posts.

## Anvi'o: a versatile software package

To get started with Anvi'o, check out this starter guide: [link][anvi-start-link]. Focus on generating the contigs database, running HMMs, profiling BAM files and merging profiles. Merging profiles will by default bin the contigs (splits) using CONCOCT. To compare multiple binning approaches, they need to be imported into Anvi'o. Follow this [guide][anvi-multi-bin-compare-link] for information on how to do that.

## Tip for importing external binning results

In the section to [import external binning results][anvi-multi-bin-compare-link] in the above tutorial, it states that the file which contains the binning results must be tab delimited, with the first column containing the contig names, and the second column containing the bin assignment. Example file:

```
contig_1	bin.1
contig_2	bin.1
contig_3	bin.2
contig_4	unbinned
contig_5	bin.2
```

In the [previous][section7-link] section, `checkm coverage` was ran, and produced a file, `coverage.tsv`. This file is a tab delimited text file that contains coverage information in each column. But more importantly, the first two columns contain the contig names and bin assignments. **So, you can subset the first two columns of this file to generate your external binning results file**.

To subset the first two columns with a tab between them, use `awk`:

```bash
awk '{printf "%s\t%s\n", $1, $2}' coverage.tsv > external_binning_results.txt
```

You may have noticed that the bin names in the example external binning results file above contains delimiters Anvi'o doesn't like (namely the `.` characters). Check to see if any of the bin names contain delimiters other than underscores - these will need to be changed to underscores. Use `sed` to achieve this:

```bash
sed 's/\./\_/g' external_binning_results.txt > tmp && mv tmp external_binning_results.txt
```

Note that different binning software will title the bin names differently, and so you should adjust the `sed` command accordingly. But in the case of MetaBAT and MaxBin2, both will delimit bin names using `.`.

You can also do this all in one line:

```bash
awk '{printf "%s\t%s\n", $1, $2}' coverage.tsv | sed 's/\./\_/g' > external_binning_results.txt
```

And now your external binning results should look like this:

```
contig_1	bin_1
contig_2	bin_1
contig_3	bin_2
contig_4	unbinned
contig_5	bin_2
```

You can now safely import this collection into Anvi'o.

## Identifying taxonomy for each bin using Kaiju

Taxonomic classification of bins is useful for determining what each bin contains, and is also useful for manual bin refinement (next section). Kaiju is a taxonomic classifier which uses Burrows-Wheeler transform to find maximum (in)exact matches at the protein level. Conveniently, there is a [guide][anvi-import-kaiju-taxa-link] on how to import taxonomic annotations from kaiju into Anvi'o. Kaiju is Linux only, and therefore must be run on the cluster. Note that for the majority of databases, Kaiju will require lots of RAM, and hence it is safest to run it on **crunch01** to avoid any problems.

To get started, a reference database must be downloaded for kaiju to search against. The full installation instructions are [here][kaiju-setup-link], and note the different databases which are available for use. However, note the size and memory requirements for each database, and choose the best one according to resources available. Note that only the files `kaiju_db.fmi` (or `kaiju_db_nr.fmi` / `kaiju_db_nr_euk.fmi`), `nodes.dmp`, and `names.dmp` are needed, and the rest can be deleted. The Mar reference databases are already downloaded to the Parfrey lab folder on the Zoology file server, found at: `/parfreylab/shared/databases/kaiju_db/mar`. It was downloaded with the following command:

```bash
# depending on the database, '-m' can be switched for other flags
makeDB.sh -m -t 12
```

If you wish to download different databases, save them to a new folder, `/parfreylab/shared/databases/kaiju_db/DB_NAME/`, where `DB_NAME` is the name of the database you are downloading. 

Note that this script may require large amounts of memory, depending on the size of the database you are downloading. Furthermore, this script doesn't work on the cluster, **and so you will need to run it on a lab machine**. After that, you can delete all the files except for the three mentioned above, and then copy those three files to `/parfreylab/shared/databases/kaiju_db/DB_NAME/`.

Continue with the [guide][anvi-import-kaiju-taxa-link] to assign and import taxonomies into Anvi'o.

## Refining bins using anvi-refine

Generally speaking, a bin with >90% completeness and <10% redundancy is considered to be a high quality draft genome, but this will not always be the case with every bin. Furthermore, high completeness and low redundancy are only two of several factors which dictate the quality of a bin, as some bins which hit these numbers may still contain erronous contigs. 

At some point, the automated binning approaches used above will need to be inspected, and manual refinement will be required for high quality bins. Manual refinement may decrease completeness; however, it's better to focus on obtaining high quality genome bins with confidence rather than focusing on hard numbers.

To get started on manual refinement, follow this [guide][anvi-meren-refine-bins-link]. This is an introduction on how to manually refine bins in Anvi'o. It focuses on lowering redundancy (contamination) by using taxonomies (which we assigned [here][kaiju-section-link]) assigned to each contig as a means to further split and refine the bins. You will be able to save and export your refined bins as well, and visualize them via `anvi-summarize`, as explained in the guide.

There is also very insightful information from this other [guide][anvi-veronica-refine-bins-link]. While the general cutoffs for a high quality draft genome is over 90% completeness and less than 10% redundancy, this guide explains the pitfalls of focusing solely on those numbers by providing examples of such bins which meet those criteria, yet still can be refined further. In addition to using taxonomy to split bins like the guide above, this guide also uses coverage, taxonomies through `BLAST`, and sequence composition based methods as additional factors to consider.

**Be sure to read through and follow along both**.

## Next step

Proceed to [section 9][section9-link].

[anvi-start-link]: http://merenlab.org/2016/06/22/anvio-tutorial-v2/
[anvi-multi-bin-compare-link]: http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-collection
[anvi-import-kaiju-taxa-link]: http://merenlab.org/2016/06/18/importing-taxonomy/#kaiju
[anvi-meren-refine-bins-link]: http://merenlab.org/2015/05/11/anvi-refine/
[anvi-veronica-refine-bins-link]: http://merenlab.org/2017/05/11/anvi-refine-by-veronika/
[kaiju-setup-link]: https://github.com/bioinformatics-centre/kaiju#creating-the-reference-database-and-index
[kaiju-section-link]: #identifying-taxonomy-for-each-bin-using-kaiju
[section7-link]: ../section_7
[section9-link]: ../section_9
