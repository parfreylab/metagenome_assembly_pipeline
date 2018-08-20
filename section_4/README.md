# Section 4: Coverage

This section obtains coverage information by mapping reads back to contigs.

## Coverage by mapping reads back to contigs

Map reads back to the assembly to get coverage information using `bowtie2`. Build the index database with `bowtie2-build` using the contigs outputted from the assembler, followed by samtools to generate a sorted and indexed `.bam` file. 

Note that multiple files can be passed into `bowtie2` by delimiting them with a comma. Alternatively, you can loop through the reads files like in the example commands below. The merged and non-merged reads can be passed through the `-U` flag.

```bash
# build the index db using bowtie2-build
bowtie2-build --threads 24 metaspades_results/contig.fasta db/spades_contigs

# bowtie2 to align reads to assembled contigs, samtools to create a .bam file. do this for EACH SAMPLE
# -1: Files with forward reads, paired with files in -2.
# -2: Files with reverse reads, paired with files in -1.
# -U: Files with unpaired reads.
for f in *_R1.fq.gz; do
    r=${f/_R1/_R2}
    m=${f/_R1/_merged}
    u=${f/_R1/_notmerged}
    out=${f/.unmapped.ecct_R1.fq.gz/}
    
    bowtie2 --threads 24 --time -x db/spades_contigs -1 $f -2 $r -U m,u | samtools view -bS -o ../$out.bam
done

# sort and index bam files
for f in *.bam; do
    out=${f/.bam/.sorted.bam}
    
    echo =======================
    echo sorting and indexing $f
    echo =======================
    
    samtools sort --threads 23 -o $out $f && samtools index $out
done
```
## Next step

Proceed to [section 5][section5-link].

[section5-link]: ../section_5