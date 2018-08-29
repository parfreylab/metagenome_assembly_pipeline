# Section 1: File preprocessing

This section deals with quality inspection, quality trimming and filtering, and optionally, removal of host sequences.

## Initial quality scores

Run `fastQC` and `multiQC` on your raw data to observe sequence quality and statistics. 

The `fastQC` output comes in 2 files, a .zip and .html file. The .html files can be viewed for each individual sample, and `multiQC` can be used to create and visualize an overall summary of quality statistics across all samples.

```bash
# -t: number of threads
fastqc -t 23 *.fastq.gz -o /path/to/output/dir && multiqc /path/to/output/dir -o /path/to/output/dir
```

Inspect the `multiQC` output. Illumina sequencing arrives a file for each lane, so this report can be used to determine if there was any problems during sequencing, or if a specific lane was problematic. If all lanes are okay, concatenate them.

```bash
# change R1 to R2 for reverse reads files
for f in *_L001_R1_*.gz; do 	# here
    l2=${f/_L001_/_L002_}
    l3=${f/_L001_/_L003_}
    l4=${f/_L001_/_L004_}
    out=${f//_L001_R1_001} 		# here
    read=R1 					# here
    
    cat $f $l2 $l3 $l4 > $f.tmp && mv $f.tmp path/to/reads/folder/${out%.fastq.gz}_$read.fastq.gz 
done
```

## Data preprocessing with BBTools

The following is mostly based off of the assembly pipeline included in the `BBMap` utils. It can be found at `path/to/bbmap/pipelines/assemblyPipeline.sh`. As a side note, the `BBTools` suite also comes with lots of documentation, which can be found under `path/to/bbmap/docs/guides/`, and `path/to/bbmap/pipelines/`. If you are using entamobea, `BBTools` will already be installed at `/Users/parfreylab/Desktop/lab_member_files/programs/bbmap`.

Remove optical duplicates, which are errors in sequences introduced by the sequencing machinery itself.

```bash
# BBMap genome assembly pipeline
for f in *_R1.fastq.gz; do
    r=${f/_R1/_R2}
    fwd=${f/.fastq/.clumpify.fq}
    rev=${r/.fastq/.clumpify.fq}
    
    echo
    echo ==============================
    echo clumpify, dedupe for sample $f
    echo ==============================
    echo
    
    clumpify.sh \
    in=$f in2=$r \
    out=../clumpify_rm_optical_dupes/$fwd \
    out2=../clumpify_rm_optical_dupes/$rev \
    dedupe \
    optical
done
```

Remove low quality reads without adding bias/salvage libraries with major positional problems in lanes. See [here](https://www.biostars.org/p/228762/) for a more detailed description. Note that this step should only be used if the dataset is not a complex metagenome, or if the samples have very low coverage.

As per author's recommendation: "However, there are some cases such as **complex metagenomes** - in which more coverage is strictly beneficial, so throwing away even low-quality reads is a bad idea.  In these cases, or any situation where very low coverage is expected, filtering will often lead to inferior results.  With high coverage, `FilterByTile` should be strictly beneficial."

```bash
# in clumpify_rm_optical_dupes folder
for f in *_R1.clumpify.fq.gz; do
    r=${f/_R1/_R2}
    fwd=${f/.clumpify.fq/.tile.fq}
    rev=${r/.clumpify.fq/.tile.fq}
    
    echo
    echo ==============================
    echo filtering by tile for sample $f
    echo ==============================
    echo
    
    filterbytile.sh \
    in=$f \
    in2=$r \
    out=../filter_by_tile/$fwd \
    out2=../filter_by_tile/$rev
done
```

Trim adapters. Note that adapters are automatically screened out by the Illumina software, so typically only a small percentage of reads are removed. For metagenomic analysis, adapters won't be functionally annotated by downstream software, but may cause spurious contigs. 

The `adapters.fa` file comes with `BBTools`. If you are using entamoeba, this can be found at `/Users/parfreylab/Desktop/lab_member_files/programs/bbmap/resources/adapters.fa`.

```bash
# in filter_by_tile folder
# Params:
# tbo (trim by overlap): also trim adapters based on paired-end overlap using BBMerge.
# tpe (trim paired ends): trim both reads to the same length.
# ordered: order of reads in = order of reads out.
# ktrim=r: trim on the 3' (right) side of reads. in this particular mode, once a reference k-mer is matched, 
# the matching k-mer and everything to the right of it will be trimmed, leaving only bases to the left.
# k=23: the k-mer length used for trimming
# mink=11: the minimum k-mer length to be considered. useful for if adapter sequences are smaller than the 
# given k, as it will use lengths from 22 down to 11 (in this case) from the right hand side of the read.
# hdist=1: minimum hamming distance (i.e. # of positions where two strings of equal length differ by).
# ftm=5 (force trim modulo): trim the right end of the read such that it's length is 0 mod 5. 
for f in *_R1.tile.fq.gz; do
    r=${f/_R1/_R2}
    out1=${f/.tile/.trimmed}
    out2=${r/.tile/.trimmed}
    
    echo
    echo ==============================
    echo adapter trimming for sample $f
    echo ==============================
    echo
    
    bbduk.sh tbo tpe ordered \
    in=$f \
    in2=$r \
    out=../bbduk_trim/$out1 \
    out2=../bbduk_trim/$out2 \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    ref=path/to/bbmap/resources/adapters.fa \
    ftm=5 
done
```

Remove low quality synthetic artifacts and spike-ins by k-mer matching. Use default parameters, as in the pipeline. Note that the reference files used (`phix174_ill.ref.fa.gz` and `sequencing_artifacts.fa.gz`) come with the installation of `BBTools`. If you are using entamoeba, these will be under `/Users/parfreylab/Desktop/lab_member_files/programs/bbmap/resources`. 

```bash
# in bbduk_trim folder
# Params:
# ref: reference fasta file for adapter sequences.
# cardinality: number of unique k-mers counted by use of the loglog algorithm 
# (approximates number of distinct elements in a multiset).
for f in *_R1.trimmed.fq.gz; do
    r=${f/_R1/_R2}
    out1=${f/.trimmed/.filtered}
    out2=${r/.trimmed/.filtered}
    
    echo
    echo ========================================================
    echo removing synthetic artifacts and spike-ins for sample $f
    echo ========================================================
    echo
    
    bbduk.sh ordered cardinality \
    in=$f \
    in2=$r \
    out=../bbduk_rm_sa_si/$out1 \
    out2=../bbduk_rm_sa_si/$out2 \
    k=31 \
    ref=path/to/bbmap/resources/phix174_ill.ref.fa.gz,path/to/bbmap/resources/sequencing_artifacts.fa.gz 
done
```

## [OPTIONAL] Removal of host sequences (decontamination)

This step is optional, as not every dataset will require this. If this does not apply, skip to [error correction][error-correction-section-link].

Host genomes can be downloaded from a variety of websites, including [EchinoBase][echinobase-link], the [NCBI][ncbi-link],

Start by building an index database with the host reference genome using `bowtie2`.

```bash
bowtie2-build --threads 23 host_scaffolds.fa bowtie2_ref_scaffolds/host_db
```

Map to the host reference genome. If reads map to the host genome, then they can be considered host sequences and subsequently removed. The unmapped reads will be the reads of interest. In the following example command, if you don't want to keep the aligned (assumed host) reads, remove the `--al-conc-gz` flag and argument.

```bash
# in bbduk_rm_sa_si folder
# Params:
# x: input database folder, constructed from the host scaffolds/contigs.
# threads: number of threads to use. 
# -1/-2: forward and reverse reads, respectively.
# al-conc-gz: paired end reads which aligned to the database (i.e. host sequences). 
# -gz appended writes gzipped files. 
# un-conc-gz: paired end reads which didn't align (i.e. host-removed sequences). 
# S: output SAM files.
db=path/to/bowtie2_ref_scaffolds/host_db
for f in *_R1_*.filtered.fq.gz; do
    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    out=${f/_R1_001.filtered.fq.gz/.mapped.fq.gz}
    outu=${out/.mapped/.unmapped}
    sam=${f/_R1_001.filtered.fq.gz/.mapped.unmapped.sam}
    
    echo
    echo ==============================
    echo mapping sample $f
    echo ==============================
    echo

    bowtie2 --very-sensitive-local \
    -x $db \
    --threads 24 \
    -1 $f \
    -2 $r \
    --al-conc-gz ../bowtie2_host_rm/$out \ # aligned .fq.gz files
    --un-conc-gz ../bowtie2_host_rm/$outu \ # unaligned .fq.gz files
    -S ../bowtie2_host_rm/sam/$sam
done
```

At this point, taxonomies of mapped (host-associated) sequences can be checked using `blast` or another read aligner to confirm the identities of each sequence. Host associated sequences can also be moved to a different location to exclude them from downstream analysis. 

## Error correction with `bbmerge.sh` in 3 phases

The explanations are mostly through the parameters, but you can always find the full documentation at `path/to/bbmap/docs/guides`.

Phase 1:

```bash
# in bowtie2_host_rm
# ecco: error correct reads that overlap rather than merging. region where two reads agree will receive 
# higher quality scores, disagree will receive lower and base call will be changed to one with higher score. 
# if scores are equal but bases are different, then "N" will be put in. 
# mix: output both merged and unmerged reads in the same file.
# vstrict: greatly reduces false positive rate, but decreases merging rate. 
# a parameter which sets a bunch of other parameters (author of bbtools suggests using strict/loose params).
# ihist: insert length histogram output file.
# overwrite: force overwrite existing files.
for f in *_R1*.fq.gz; do
    r=${f/_R1/_R2}
    out1=${f/_R1/.ecco_R1}
    out2=${r/_R2/.ecco_R2}
    hist=${f/.fq.gz/.hist.txt}
    
    echo 
    echo =============================== 
    echo error correction 1 for sample $f 
    echo =============================== 
    echo 
    
    bbmerge.sh ecco mix vstrict ordered \
    in=$f \
    in2=$r \
    out=../bbmerge_ecco_1/$out1 \
    out2=../bbmerge_ecco_1/$out2 \
    ihist=../bbmerge_ecco_1/hists/$hist \
    overwrite=true
done
```

Phase 2:

```bash
# in bbmerge_ecco_1
# ecc: error correct reads, requiring multiple passes.
# passes=4: number of error correction passes (4).
for f in *_R1*.fq.gz; do
    r=${f/_R1/_R2}
    out1=${f/.ecco_R1/.eccc_R1}
    out2=${r/.ecco_R2/.eccc_R2}
    
    echo 
    echo =============================== 
    echo error correction 2 for sample $f 
    echo =============================== 
    echo 
    
    clumpify.sh ecc passes=4 reorder \
    in=$f \
    in2=$r \
    out=../bbmerge_eccc_2/$out1 \
    out2=../bbmerge_eccc_2/$out2 \
    overwrite=true
done
```

Phase 3: Error correction and extension of reads. Use ~1/3 of the read length as the k-mer length for error correction.

```bash
# in bbmerge_eccc_2
for f in *_R1*.fq.gz; do
    r=${f/_R1/_R2}
    out1=${f/.eccc_R1/.ecct_R1}
    out2=${r/.eccc_R2/.ecct_R2}
    
    echo 
    echo =============================== 
    echo error correction 3 for sample $f 
    echo =============================== 
    echo 
    
    tadpole.sh ecc k=31 ordered \
    in=$f \
    in2=$r \
    out=../bbmerge_ecct_3/$out1 \
    out2=../bbmerge_ecct_3/$out2 \
    overwrite=true
done
```

## Merge overlapping reads, quality filter unmerged ones

```bash
# in bbmerge_ecct_3
# strict: decrease false positive rate and merging rate.
# rem (requireextensionmatch): do not merge if predicted insert size differs before and after extension. 
# will attempt to merge raw reads, then extend and try to merge again. If both merges give same result, 
# then reads will be merged. If they give different results, the reads will only be merged if the raw reads 
# gave no solution but the extended reads gave a solution indicating the raw reads don't overlap, or if the 
# raw reads DID have a solution and no extension was possible. 
# extend2: Extend reads this much only after a failed merge attempt.
for f in *_R1*.fq.gz; do
    r=${f/_R1/_R2}
    out=${f/.ecct_R1/.merged}
    outn=${f/.ecct_R1/.notmerged}
    hist=${f/.ecct_R1.fq.gz/.hist.txt}
    
    echo 
    echo =============================== 
    echo merging paired-end sample $f 
    echo =============================== 
    echo 
    
    bbmerge-auto.sh strict rem ordered \
    in=$f \
    in2=$r \
    out=../bbmerge_merged_k31/$out \
    outu=../bbmerge_merged_k31/$outn \
    k=31 \
    extend2=30 \
    ihist=../bbmerge_merged_k31/hist/$hist
done
```

Generate k-mer depth distributions.

```bash
# in bbmerge_merged_k31
for f in *.merged.fq.gz; do
    o=${f/.merged.fq.gz/.khist.txt}
    o2=${f/.merged.fq.gz/.peaks.txt}
    
    khist.sh in=$f khist=khist/$o peaks=khist/$o2
done
```

Quality trim the unmerged reads. Here, the trim quality (`trimq`) and minimum length (`minlen`) cutoffs can be adjusted to fit your data.

```bash
# qtrim: Trim read ends to remove bases with quality below trimq.
# Performed AFTER looking for kmers. Values: 
# rl (trim both ends), 
# f (neither end), 
# r (right end only), 
# l (left end only),
# w (sliding window).
# trimq: Regions with average quality BELOW this will be trimmed, if qtrim is set to something other than f. 
# Can be a floating-point number like 7.3.
for f in *.notmerged.fq.gz; do
    out=${f/.notmerged/.notmerged.trimmed}
    
    echo 
    echo =============================== 
    echo trimming non-merged sample $f 
    echo =============================== 
    echo 
    
    bbduk.sh in=$f out=../bbduk_nm_qtrim20_ml70/$out qtrim=r trimq=20 minlen=70 ordered
done
```

## Merging samples for assembly (IDBA-UD)

In the next section, these preprocessed reads will be assembled into contigs and scaffolds using `metaSPAdes` and `IDBA-UD`. `metaSPAdes` accepts multiple samples as input, but `IDBA-UD` does not. Hence, we need to concatenate and create the following files as input into `IDBA-UD`:

1. concatenated merged reads
2. concatenated non-merged forward and reverse quality filtered reads
3. concatenated forward reads, prior to merging 
4. concatenated reverse reads, prior to merging  

Files 2 and 3 can be made from after phase 3 error correction (just before merging). Follow a similar command below for each file.

```bash
# in bbmerge_merged_k31 folder 
for f in *.merged.fq.gz; do
    cat $f >> assembly_data/all_data_merged.fq.gz
done

# in bbmerge_merged_k31 folder
for f in *.notmerged.fq.gz; do
    cat $f >> assembly_data/all_data_notmerged.fq.gz
done

# in bbmerge_ecct_3 folder
for f in *.ecct_R1.fq.gz; do
	cat $f >> assembly_data/all_data_forward_reads.fq.gz
done

# in bbmerge_ecct_3 folder
for f in *.ecct_R2.fq.gz; do
	cat $f >> assembly_data/all_data_reverse_reads.fq.gz
done
```

## Next step

Proceed to [section 2][section2-link].

[section2-link]: ../section_2
[error-correction-section-link]: #error-correction-with-bbmergesh-in-3-phases
[ncbi-link]: https://www.ncbi.nlm.nih.gov/
[echinobase-link]: http://www.echinobase.org/Echinobase/