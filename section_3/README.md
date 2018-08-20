# Assembly evaluation

This step evaluates the quality of the assemblies from `metaSPAdes` and `IDBA-UD` from the previous step. 

## Assembly evaluation with metaQUAST

Use `metaQUAST` evaluate basic statistics and overall quality of the assembly. Such statistics include N50, L50, total assembly size, etc. Generally, a higher N50 value and a lower L50 value implies a more contiguous assembly, which is better than a fragmented one. However, size statistics are only one of several criteria you should use to evaluate the overall quality of an assembly. A high N50 may not mean much if all the contigs are spurious, or are chimeric. 

Gene finding can also be specified by passing the `-f` flag. With `metaquast.py`, `MetaGeneMark` will be used as the gene finder. If no reference genomes are provided, `metaQUAST` will blast the contigs against the `SILVA` 16s database, and pull out similar reference genomes past a certain threshold.

```bash
metaquast.py --threads 16 -f -o quast_spades/ path/to/metaspades/output/contigs.fasta 
metaquast.py --threads 16 -f -o quast_idba/ path/to/idba/out/contigs.fa
```

Compare the two assemblies. Which assembly has the better size statistics? Which one has more longer contigs? Which one had more fully unaligned contigs? These questions may help you in evaluating which assembly is better. Alternatively, the pipeline can be continued with both assemblies and the results can be compared.

## Next step

**NOTE**: the following sections will assume you chose the `metaSPAdes` assembly and associated files (for no reason other than to simplify the commands). Proceed to [section 4][section4-link].

[section4-link]: ../section_4