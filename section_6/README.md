# Binning

In this section, the contigs obtained from the previous two sections will be binned using multiple binning programs. Each method listed here can be used together, or a single one can be used; however, it's recommended to try each one to see which one gives the best results. In downstream analysis (in An'vio), multiple binning approaches can be compared to see how they compare.

## MetaBAT

Binning with `metabat2`. MetaBAT uses tetranucleotide frequencies (TNF) and abundance probabilities. TNFs of different genomes are discriminated against using likelihood of inter + intra species Euclidian distance from 1414 microbial genomes. This distance is known as the tetranucleotide frequency distance probability (TDP).

For each pair of contigs in a metagenome assembly, MetaBAT calculates their probabilistic distances based on tetranucleotide frequency (TNF) and abundance (i.e., mean base coverage), then the two distances are integrated into one composite distance. All the pairwise distances form a matrix, which then is supplied to a modified k-medoid clustering algorithm to bin contigs iteratively and exhaustively into genome bins.

Do the following step using the ZCU, since the `runMetaBat.sh` script requires usage of Linux binaries. MetaBAT requires a sorted and indexed BAM file (obtained from [mapping reads back to contigs][section4-link]) in order to generate a depth file.

The authors of MetaBAT also have a short [tutorial][metabat2-binning-tutorial-link] on how to refine parameters passed to `metabat2`. Note that this requires use of `CheckM`, `R` and some custom scripts.

```bash
# -t: number of threads to use
# -v: verbose mode (good for seeing progress)
# --unbinned: keep all unbinned sequences in a separate fasta file. may be useful later on
# -m: min contig length (default 2500)
runMetaBat.sh -t 8 -v --unbinned -m 2000 filtered_contigs.fasta sorted.bam
# runMetaBat.sh doesn't let you specify the output folder name, so rename it like so
mv filtered_contigs.fasta.metabat-bins2000 metabat_bins
```

## MaxBin2

Binning with MaxBin2. MaxBin2 is an alternative to MetaBAT, and uses tetranucleotide frequencies and scaffold/contig coverages. Furthermore, it uses single copy marker genes and an Expectation-Maximization (EM) algorithm to populate the N genomic bins. Essentially, the EM algorithm calculates a probability for each scaffold/contig to belong to each genomic bin. Scaffolds/contigs are then assigned to the bin with the highest probability. MaxBin2 will also recursively split bins based on the number of marker genes it finds in the bin.

Do the following step on the ZCU, as MaxBin2 has been tested and written under Linux platforms. It may be possible to use on macOS, but there are no guarantees (installation is also harder, as the automation script won't work on macs). Depending on whether or not a tab-delimited abundance file is available, MaxBin2 can take this abundance file or reads can be input for mapping with Bowtie2. See the [README][maxbin2-readme-link] for more details. 

```bash
# with the abundance list
perl run_MaxBin.pl -contig contigs.fa -abund tab-delimited-abundances.txt -out path/to/out/out
```

If there are multiple samples, then a reads list can be provided instead. The list is a separate file with each line as a path to a single reads file.

```bash
# with reads. this is the simplest option if running for the first time
perl run_MaxBin.pl -contig contigs.fa -reads_list reads_list.txt
```

## CONCOCT via An'vio

CONCOCT is yet another binning tool, but is shipped with An'vio (later used for manual bin refinement, visualization, and summaries). CONCOCT uses sequence composition and coverage across multiple samples for binning. The Meren lab (developers of An'vio) have great resources on their website on how to use An'vio and all of its functionality, so any steps regarding An'vio in this tutorial will be redirected to their [website][meren-lab-anvio-link]. Check the section [Binning, taxonomy and bin befinement in An'vio][section8-link] for more information.

## Next step

Proceed to [section 7][section7-link].

[section4-link]: ../section_4
[metabat2-binning-tutorial-link]: https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices
[maxbin2-readme-link]: https://downloads.jbei.org/data/microbial_communities/MaxBin/README.txt
[section7-link]: ../section_7
[meren-lab-anvio-link]: http://merenlab.org/software/anvio/
[section8-link]: ../section_8
