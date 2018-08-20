# Section 0: Prelude

The overall process and technical challenges you may face during analysis of metagenomic data can be hard at times. Some tools have great documentation, are up-to-date, and
have good support, while others have less documentation and are less usable. There are also tons of tutorials focusing on the overall process of metagenome analysis - but with all
these different tools and documentations, it can get overwhelming. To try and make things simpler, here is a list of guides I referred to the most while writing this pipeline:

- [The simple fool's guide to population genomics via RNA-Seq: an introduction to high-throughput sequencing data analysis.][stanford-overall-workflow-link] A good guide which focuses on the overall process of data analysis, with in-depth explanations of each step.  
- [An introduction to genome assemblies.][intro-genome-assemblies-link] Another guide to look at to understand the overall workflow of a genome assembly pipeline.
- [Metagenomics wiki.][metagenomics-wiki-link] A page which compiles definitions of frequently used terms, tools, and some workflows for specific steps of a pipeline.
- [JGI data preprocessing.][jgi-data-preprocess-link] How the Joint Genome Institute (JGI) prepares their data for analysis.
- [Biostars][biostars-link] A forum for help with bioinformatics.
- [SEQanswers][seq-answers-link] Another forum for help with bioinformatics.

## Adjusting parameters

For questions like: What parameters should I adjust? What k-mer size should I use? Am I being too stringent with filtering? Often times there is no one universal answer for best practices, as each dataset is different. But generally, the tools used in this pipeline have default values for each parameter. The authors of each tool likely had good reasons to arrive at the values they did, and so it's best not to stray too far from those values unless you have a compelling reason to. It's always a good idea to use Google search, Biostars, and SEQanswers (linked above) as well to get more opinions if you are unsure. 

## Obtaining the data

Most data will be found on the Parfrey lab folder under the Zoology file server. If you are working with the data on a local machine, use `rsync` to copy it locally. Here is an example command with an [explanation][rsync-explanation-link] of the flags:

```bash
rsync -chavzP --stats [source] [target]
```

Which may look more like this for our use case:

```bash
rsync -chavzP --stats your.login@zoology.ubc.ca:/path/to/data /path/to/data/folder/locally
```

## Next step

Proceed to [section 1][section1-link].

[stanford-overall-workflow-link]: http://sfg.stanford.edu/guide.html
[intro-genome-assemblies-link]: https://sschmeier.github.io/bioinf-workshop/genome-assembly/doc/GenomeAssembly_sschmeier.pdf
[metagenomics-wiki-link]: http://www.metagenomics.wiki/pdf/definition
[jgi-data-preprocess-link]: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/data-preprocessing/
[biostars-link]: https://www.biostars.org/
[seq-answers-link]: http://seqanswers.com/forums/
[section1-link]: ../section_1/
[rsync-explanation-link]: https://explainshell.com/explain?cmd=rsync+--chavzP+--stats+%5Bsource%5D+%5Btarget%5D