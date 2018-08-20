# Section 5: Filtering contigs

This step filters out any unwanted contigs from the assembly by aligning them to a database (`NT` in this example). **NOTE** that this step is **OPTIONAL**, and can be [skipped][section6-link] if it does not apply.

## [OPTIONAL] Filtering out unwanted contigs

Depending on your data, you may want to filter out contigs which are not of interest. Run this on the cluster, as `blastn` is installed globally and there is a copy of `NT` in the Parfrey lab file server folder. There may be other databases you would like to `BLAST` against - however, you will need to install these yourself.

NOTE: this command must be run in the `bash` shell. To do this, type `bash` in the terminal when logged onto the ZCU. Then, set the `$BLASTDB` variable by adding it to `.bashrc`. See [this][set-up-blast-db-link] blog post and [this][biostars-blast-link] thread for more installation setup. The first link has information about what the headers mean of each column in the tabular `BLAST` output file, and the second link explains how to include taxonomies for each hit.

```bash
# add this to your .bashrc
export BLASTDB="/parfreylab/shared/databases/NCBI_NT/NT_06.08.2018"
```

Run `blastn` like the following example. Note that the output format should keep the same columns and order as specified in the example command so `MEGAN` will be able to parse it. 

Depending on how distantly related your contigs are to the `NT` reference, you should consider changing the percent identity (`-perc_identity`) and e-value (`-evalue`) cutoffs. Also, consider changing the flag `-num_threads` depending on the number of jobs a machine has running.

```bash
blastn \
-db nt \
-query /parfreylab/kevchan/assembly/idba_ud_contigs/contig.cec.fasta \
-num_threads 8 \
-max_hsps 1 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames' \
-out blast_filter_contigs/cec_idba_contigs_hit.txt \
-perc_identity 90 \
-max_target_seqs 10 \
-evalue 0.005
```

Inspect the resulting file to confirm the results are what is expected. Otherwise, some of the above parameters may need to be changed. If the results look good, then use `MEGAN` to assign taxonomies to the sequences. `MEGAN` uses a lowest common ancestor algorithm to assign a taxonomy to each sequence, and can be viewed through their GUI. Import the blast tabular output file and contigs so that they can be extracted. Remove any unwanted contigs.

Once you are in the `MEGAN` interface, import the `BLAST` results by going to `File > Import from BLAST` (or `SHIFT-CMD-N`, the shortcut). Import the contigs files as well in the same interface under long reads mode, so that once `MEGAN` computes the LCA for each contig, you can select which contigs to keep for downstream analysis. 

If you don't have one yet, down the latest copy of the accession mapping file, titled `nucl_acc2tax-MONTHYEAR.abin.zip` from the `MEGAN` [website][megan-website-link]. Unzip the file, and return to the `MEGAN` interface. Go to the `Taxonomy` tab, click on 'Load Accession mapping file', and input the newly unzipped `nucl_acc2tax` file. This will allow `MEGAN` to assign the contigs a taxonomy. Leave the rest as defaults, and then click `Apply`. 

In the resulting view, the nodes can be clicked from within the GUI. The level of classification can also be adjust by clicking on the `Rank...` tab on the top pane. Once you have found your desired contigs, they can be exported by right clicking on a node, and clicking `Extract Reads...`.

## Next step

Proceed to [section 6][section6-link].

[section6-link]: ../section_6
[set-up-blast-db-link]: https://iamphioxus.org/2018/01/08/local-installation-of-ncbi-blast-together-with-the-nr-and-taxonomy-database
[biostars-blast-link]: https://www.biostars.org/p/76551/
[megan-website-link]: http://ab.inf.uni-tuebingen.de/data/software/megan6/download/welcome.html
