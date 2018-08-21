# Section 10: Functional annotation with PROKKA and Krona plots

This section functionally annotates bins from the previous section with PROKKA, and generates Krona (interactive pie chart) plots.

## Functional annotation with PROKKA and Krona plots

This section of the pipeline closely follows a tutorial from a metagenomics workshop held in Uppsala in 2014. We will focus on the section "Functional Annotation", which is found [here](http://metagenomics-workshop.readthedocs.io/en/latest/annotation/functional_annotation.html). Since this tutorial is from 2014, the reference database files used to convert from Enzyme Commission (EC) numbers to hierarchies/categories to pathways are likely not up to date. Nonetheless, this should give you a good visualization of your data.

[PROKKA](https://github.com/tseemann/prokka) is a rapid and powerful whole genome annotation tool, which serves as pipeline to generate a variety of files. PROKKA can be used to generate [Krona](https://github.com/marbl/Krona/wiki) plots with annotations using 3 databases: COG, KEGG, and Metacyc.

A couple things to note:

- Skip the section "Taxonomic annotation" as this has been done with kaiju in the [previous section][kaiju-annotation-link].

- Under the section "Mapping reads and quantifying genes", work with the `.bam` files created in the previous section of this guide, [Coverage by mapping reads back to contigs][section4-link].

- The bash command to extract COG identifiers doesn't work with newer versions of PROKKA. Instead, use the two commands below. Replace `PROKKA_FILE` with the actual filename(s) instead. If there are multiple `.gff` files, run these commands in a loop.

```bash
# lines with no EC number. COG id will be in third field
egrep "COG[0-9]{4}" PROKKA_FILE.gff | \
    cut -f9 | grep -v "eC_number=" | \
    cut -f 1,3 -d ';' | \
    gsed 's/ID=//g' | \
    gsed 's/;dbxref=COG:/\t/g' | \
    grep COG > cog/PROKKA_FILE.cog

# lines WITH EC number. COG id in fourth field
egrep "COG[0-9]{4}" PROKKA_FILE.gff | \
    cut -f9 | grep "eC_number=" | \
    cut -f 1,4 -d ';' | \
    gsed 's/ID=//g' | \
    gsed 's/;dbxref=COG:/\t/g' | \
    grep COG >> cog/PROKKA_FILE.cog
```

- After removing duplicates with Picard `MarkDuplicates`, the unpaired reads need to be removed as well, otherwise `htseq-count` will not work with a `.bam` file that contains both paired and unpaired reads. Use `samtools view -bf 1` as in the following example:

```bash
for f in *.bam; do 
    name=$(echo $f | cut -d '.' -f 1)
    echo processing $name 
    samtools view -bf 1 $f | htseq-count -r pos -t CDS -f bam - ../$name/PROKKA.$name.gtf > $name.count 
done
```

## Next step

None! You've made it to the end. Congratulations!

[kaiju-annotation-link]: ../section_8#identifying-taxonomy-for-each-bin-using-kaiju
[section4-link]: ../section_4
