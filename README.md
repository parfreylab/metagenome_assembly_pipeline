# A pipeline for metagenomic assembly

A pipeline to go from raw sequencing data to high quality bins and pretty plots.

## Introduction

This is a simple pipeline for analyzing metagenomic data, starting from raw reads to assembly, to binning, to functional analysis. Throughout the pipeline, there 
will be example commands at each step - however; replace the generic filenames and paths with your own data. 

Prior to getting started, look over the [software list](#exhaustive-list-of-software-used) and ensure all of it is installed. If you are using the lab machine 'entamoeba', the majority of software will be installed locally. For steps which require the ZCU (Zoology computing unit), you will need to install the software to your own folder if it is not already installed globally.

You can preview each step of the pipeline under the [overview](#overview) section.

## Exhaustive list of software used

Note that many of these tools can be installed through conda or homebrew. Try those options first before a manual install.

### Quality metrics/visualization 

- [FastQC][fastqc-link]
- [MultiQC][multiqc-link]

### Data preprocessing

- [BBTools][bbtools-link]

### Read alignment (host removal, coverage, gene quantification)

- [bowtie2][bowtie2-link]
- [Picard][picard-link] 
- [htseq][htseq-link] 

### Assemblers

- [IDBA-UD][idba-ud-link] 
- [MetaSPAdes][metaspades-link]<sup>1,2</sup> (**NOTE:** may or may not work with your data, fix should be rolled out in v3.12.1) 

### Assembly evaluation

- [MetaQUAST][metaquast-link]

### Binning

- [CONCOCT][concoct-link] (built-in in An'vio) 
- [MetaBAT][metabat-link]<sup>1</sup> 
- [MaxBin2][maxbin2-link]<sup>1</sup> 

### Binning quality evaluation

- [CheckM][checkm-link]<sup>1</sup> 
- [An'vio][anvio-link] 

### Taxonomy assignment

- [Kaiju][kaiju-link]

### Manual binning refinement

- [An'vio][anvio-link] 

### Functional annotation

- [PROKKA][prokka-link] 
- [MinPath][minpath-link] 
- [MetaPathways][metapathways-link] 
- [GhostKOALA][ghostkoala-link] 

### Functional annotation databases

- [COGs][cogs-link] 
- [KEGG][kegg-link] 

### Other

- [Samtools][samtools-link] 
- [Krona][krona-link] 
- [blast][blast-link] 
- [MEGAN][megan-link] 

### Legend

<sup>1</sup>: Linux/cluster use only <br/>
<sup>2</sup>: High memory usage

## Overview

## Issues

Any issues or clarification needed can be raised by creating a new issue in this repository. Alternatively, you can email me at kevchan1@alumni.ubc.ca.

[fastqc-link]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[multiqc-link]: http://multiqc.info/
[bbtools-link]: https://sourceforge.net/projects/bbmap/
[bowtie2-link]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[picard-link]: https://broadinstitute.github.io/picard/
[htseq-link]: https://htseq.readthedocs.io/en/release_0.10.0/overview.html
[idba-ud-link]: https://github.com/loneknightpy/idba
[metaspades-link]: http://cab.spbu.ru/files/release3.12.0/manual.html
[metaquast-link]: http://quast.bioinf.spbau.ru/manual.html
[concoct-link]: https://concoct.readthedocs.io/en/latest/
[metabat-link]: https://bitbucket.org/berkeleylab/metabat
[maxbin2-link]: https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html
[checkm-link]: https://github.com/Ecogenomics/CheckM/wiki
[anvio-link]: http://merenlab.org/software/anvio/
[kaiju-link]: https://github.com/bioinformatics-centre/kaiju
[prokka-link]: https://github.com/tseemann/prokka
[minpath-link]: http://omics.informatics.indiana.edu/MinPath/
[metapathways-link]: https://github.com/hallamlab/metapathways2
[ghostkoala-link]: https://www.kegg.jp/ghostkoala/
[cogs-link]: https://www.ncbi.nlm.nih.gov/COG/
[kegg-link]: https://www.kegg.jp/kegg/
[samtools-link]: http://www.htslib.org/
[krona-link]: https://github.com/marbl/Krona/wiki
[blast-link]: https://www.ncbi.nlm.nih.gov/books/NBK279684/
[megan-link]: https://ab.inf.uni-tuebingen.de/software/megan6

