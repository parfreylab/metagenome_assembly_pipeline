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

## Running programs on the cluster using screen

TL;DR - try to use a screen wherever you can.

Whenever you use the cluster, whether it be for mandatory steps in this pipeline or for your own convenience, it is strongly recommended to run your programs in a *screen*. In short, a screen allows you to run a program even if your terminal session is not active. For most programs you encounter (and all the ones used in this pipeline), you will need an active connection to the program (i.e. your machine must be on and running). This tends to be problematic for steps which can take a long time - take the [assembly][section2-link] section for example: I've experienced assembly with metaSPAdes taking over a week to finish. That's a long time for your machine and for you to be connected to the cluster continuously. If you run metaSPAdes in a screen however, you can "detach" from the screen, close your machine, and metaSPAdes will continue to run on the cluster until completion (assuming no errors occurred). Here are some basic commands to work in a screen:

This command starts a new screen with session name `YOUR_SESSION_NAME`: 

```bash
screen -S YOUR_SESSION_NAME
```

To detach from a screen **and not kill the program you are running**, press `CTRL + D`. To get a list of all screens on the current machine, type:

```bash
screen -ls
```

A "detached" screen means you are not currently connected to it, while an "attached" screen means you are currently "inside" that screen. 

To reattach to a screen, do the following:

```bash
# if there is a single screen
screen -r

# if there are multiple screens
screen -r YOUR_SESSION_NAME
```

Note that multiple screens can be started on the same computer, so you may need to pass in the session name. Also note that you don't need to type in the full session name, rather just enough letters for the session name to be unique. 

To terminate a screen, kill the running program (if there is one) with `CTRL + C` or `CTRL + \`. Then type in `exit`, and you should see something like `[screen is terminating]`. After that, the screen will be terminated.

## Next step

Proceed to [section 1][section1-link].

[stanford-overall-workflow-link]: http://sfg.stanford.edu/guide.html
[intro-genome-assemblies-link]: https://sschmeier.github.io/bioinf-workshop/genome-assembly/doc/GenomeAssembly_sschmeier.pdf
[metagenomics-wiki-link]: http://www.metagenomics.wiki/pdf/definition
[jgi-data-preprocess-link]: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/data-preprocessing/
[biostars-link]: https://www.biostars.org/
[seq-answers-link]: http://seqanswers.com/forums/
[section1-link]: ../section_1/
[section2-link]: ../section_2/
[rsync-explanation-link]: https://explainshell.com/explain?cmd=rsync+--chavzP+--stats+%5Bsource%5D+%5Btarget%5D