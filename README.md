# uORF Annotator v. 0.4
*uORF Annotator* is the tool for annotation of upstream translation initiation sites in upstream open reading frames (uORFs) of human genome, which predicted by [uBert model](https://github.com/skoblov-lab/uBERTa).

New in v. 0.4:
* all effects are now annotated with respect to the uORF-specific isoform
* major issues in BED generation and main CDS effect annotation fixed

## Conda environment
Install all dependencies from `requirements.yml` as new conda environment.
```
conda env create -f requirements.yml
```
## Required input data
* VCF file of variants for further annotation
* BED file with available uORFs (`sorted.v3.bed` in this repository)
* GTF file with genomic features annotation

It is highly recommended to run the tool **with a GTF file containing uORF-matching transcript isoforms for each gene** (`combined_uorf.gtf` in this repository). If a GTF file lacks matching transcript IDs, main CDS effect will not be properly annotated.
## Run example
```
python uORF_annotator.py \
    -i <input_variants.vcf> \
    -b <input_uORFs.bed> \
    -g <input_annotation.gtf> \
    -f <human_genome.fasta> \
    -ot <output.tsv> \
    -ov <output.vsc> \
    -ob <output.bed> \
    -atg_only
```
## Output formats specification
### tab-separated (tsv) file
Each row represents annotation of a single variant in particular uORF (per uORF annotation).
#### Fields in the TSV output
1) #CHROM - contig name  
2) POS - position  
3) REF - reference allele
4) ALT - alternative allele
5) INFO - old INFO field from inputed `.vcf` file  
6) orf_start - uORF start position
7) orf_end - uORF end position
8) bed_anno - annotation from 4th field of input `.bed` file
9) strand - strand dirrection defined as + (forward) or - (reverse).
10) n_exons - number of exones in uORF  
11) exons_sizes - number of exones in uORF  
12) exons_starts - sequence of exons start positions  
13) utid - transcript id  
14) overlapping_type - type of uORF (overlapping/non-overlapping) from input `.bed` file  
15) codon_type - ATG or non-ATG start codon in uORF 
16) dist_from_orf_to_snp - distance from start of uORF to variant position  
17) technical column
18) technical column
19) symbol - schematic indication of the variant event  
20) consequence - type of effect 
21) main_cds_effect - effect of uORF variant on the main CDS (annotated for `frameshift` and `stop lost` variants only)
22) additional_info - comment field
23) in_known_ORF - a binary flag indicating if a variant affects main ORF of the gene
24) in_known_CDS - a binary flag indicating if a variant also affects main CDS of the gene
### The Variant Call Format (vcf)
Add uBERT field to INFO fields of input vcf file.
#### FORMAT:
ORF_START|ORF_END|ORF_SYMB|ORF_CONSEQ|overlapped_CDS|utid|overlapping_type|codon_type
