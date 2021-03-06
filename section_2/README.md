# Section 2: Assembly

This section deals with assembling the preprocessed reads using `metaSPAdes` and `IDBA-UD`.

## Assembly with metaSPAdes

Run `metaSPAdes` on the cluster. You **MUST** use `crunch01` since `metaSPAdes` is resource intensive, specifically in RAM. It won't work locally, even on entamoeba	 (trust me, I've tried).

Specify the forward, reverse, merged, and non-merged `fq.gz` files either through the command line, or by using a YAML dataset. Note however, that only up to 9 paired-end/mate-pair libraries can be specified via the command line. If more needs to be specified, use a `YAML` dataset (example shown below).

```bash
# metaspades input files using the command line
# -1: forward reads
# -2: reverse reads
# --merged: merged reads
# -s: single reads (quality filtered non-merged)
metaspades.py \
-1 path/to/your/forward_reads.fq.gz \
-2 path/to/your/reverse_reads.fq.gz \
--merged path/to/your/merged_reads.fq.gz \
-s path/to/your/notmerged.fq.gz \
--threads 3 \
-o path/to/metaspades_output/
```

A `YAML` dataset file can also be passed in to avoid using multiple flags. Example file:

```yaml
[
  {
    orientation: "fr",
    type: "paired-end",
    left reads: [
      "path/to/your/forward_reads.fq.gz"
    ],
    right reads: [
      "path/to/your/reverse_reads.fq.gz" 
    ],
    merged reads: [
      "path/to/your/merged_reads.fq.gz"
    ],
    single reads: [
      "path/to/your/notmerged.fq.gz"
    ]
  }
]
```

```bash
# metaspades called with YAML file
metaspades.py \
-t 4 \
--dataset path/to/metaspades_cluster_dataset.yaml \
-o . 
```

## Assembly with IDBA-UD 

For assembly with `IDBA-UD`, use the script `fq2fa` which comes with `IDBA-UD` to interleave the paired end `fastq` files, and convert them to `FASTA` format. `IDBA-UD` also uses uncompressed files only, so you must decompress any `.gz` files first before passing them to the script. Make sure there is enough disk space available before decompressing the files!

```bash
# interleaving forward and reverse reads
# --paired: if the reads are paired-end in two files
# --filter: filter out reads containing 'N'
fq2fa --merge --filter \
path/to/your/assembly_data/all_data_forward_reads.fq \
path/to/your/assembly_data/all_data_reverse_reads.fq \
paired_reads.fa
```

To pass merge and non-merged quality filter reads, concatenate the two files together, then run `fq2fa` on it. Pass in this file to the `-l` flag.

```bash
# concatenate the files
cat path/to/your/assembly_data/all_data_merged.fq.gz path/to/your/assembly_data/all_data_notmerged.fq.gz > merged_and_not_merged_reads.fq.gz
# decompress (warning: may be a large file)
gunzip merged_and_not_merged_reads.fq.gz
# 'interleaving' merged and non-merged reads
fq2fa --filter path/to/your/assembly_data/merged_and_not_merged_reads.fq merged_and_not_merged_reads.fa
```

Run `IDBA-UD` like so:

```bash
# -r: fasta read file (<= 600)
# -l: fasta long read files (> 600)
# --num_threads: number of threads to use
idba_ud \
-r paired_reads.fa \
-l merged_and_not_merged_reads.fa \
--num_threads 23 \
-o path/to/output
```

## A word on merging assemblies

There are programs out there which take sets of assemblies, and attempts to merge the assemblies to produce a "superassembly", which theoretically should result in a higher quality, more contiguous (by making supercontigs) meta-assembly. While this step would be ideal, the results may only add value to your analysis if your assemblies were high quality in the first place, and may not consistently improve quality either ([source][merging-assemblies-benchmark-link]). However, it is your analysis, and so you are welcome to try this step to see if it improves your assemblies. 

As a side note, the paper compares seven meta-assembly programs: CISA, GAA, GAM_NGS, GARM, Metassembler, MIX, and ZORRO. I (the author, Kevin) have tried and failed at installing [GARM][garm-link] and [Metassembler][metassembler-link] both on lab machines and on the Zoology cluster. Both appear to have dependencies which don't work on macOS and with the version of Linux on the cluster, so use at your own discretion.

As a final note, the last updates for GARM and Metassembler came out in 2015, so it appears that the authors may not be actively supporting the software.

## Next step

Proceed to [section 3][section3-link].

[section3-link]: ../section_3
[merging-assemblies-benchmark-link]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1213-3
[garm-link]: http://garm-meta-assem.sourceforge.net/
[metassembler-link]: http://garm-meta-assem.sourceforge.net/
