# Section 4: Coverage

This section obtains coverage information by mapping reads back to contigs.

## Preparing the FASTA file

Map reads back to the assembly to get coverage information using `bowtie2`. Build the index database with `bowtie2-build` using the contigs outputted from the assembler, followed by samtools to generate a sorted and indexed `.bam` file.

Before building the database, simplify the deflines (headers) of the fasta input file. This is necessary, as An'vio, used in a [later step][section8-link], requires simple deflines in the FASTA file to import. The names in the `.bam` file must also match, hence we simplify first before mapping to avoid repeating this step. Read the section on file formatting in the An'vio tutorial [here][anvi-starter-tutorial-link] for more information. To simplify the headers, run the script `anvi-script-reformat-fasta`. The `-l` flag specifies the minimum length of a contig to keep. By passing in `0`, we are keeping all the contigs. `--simplify-names`, as expected, simplifies the contig names:

```bash
anvi-script-reformat-fasta contigs.fasta -o contigs-fixed.fasta -l 0 --simplify-names 
```

At this point you can overwrite `contigs.fasta` with `contigs-fixed.fasta`:

```bash
mv contigs-fixed.fasta contigs.fasta
```

However, you may prefer to keep the original `contigs.fasta`, which is fine as well. The rest of this section assumes you overwrote `contigs.fasta`.

## Coverage by mapping reads back to contigs

Note that multiple files can be passed into `bowtie2` by delimiting them with a comma. Alternatively, you can loop through the reads files like in the example commands below. The merged and non-merged reads can be passed through the `-U` flag.

```bash
# build the index db using bowtie2-build
bowtie2-build --threads 24 metaspades_results/contigs.fasta db/spades_contigs

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

[anvi-starter-tutorial-link]: http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file
[section8-link]: ../section_8
[section5-link]: ../section_5