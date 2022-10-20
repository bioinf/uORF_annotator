# uORF Annotator v. 0.7
*uORF Annotator* is the tool for annotation of the functional impact of genetic variants in upstream open reading frames (uORFs) in the human genome, which are predicted by [uBert model](https://github.com/skoblov-lab/uBERTa).

New in v. 0.7:
* annotation of `stop_gained` variants with potentional activating effect on main ORF (`overlap_removed`)
* generation of a separate VCF file as a companion for bED visualization
* re-organization of TSV output and VCF fields
* gnomAD constraint metrics annotation added under `-gc` option. 

## Conda environment
Install all dependencies from `requirements.yml` as new conda environment.
```
conda env create -f requirements.yml
```
## Required input data
* VCF file of variants for further annotation
* BED file with available uORFs (`sorted.v3.bed` in this repository)
* GTF file with genomic features annotation
* \[optional\] TSV file with gene-level gnomAD constraint statistics

It is highly recommended to run the tool **with a GTF file containing uORF-matching transcript isoforms for each gene** (`combined_uorf.v4.gtf` in this repository). If a GTF file lacks matching transcript IDs, main CDS effect will not be properly annotated.
## Run example
```
python uORF_annotator.py \
    -i <input_variants.vcf> \
    -b <input_uORFs.bed> \
    -g <input_annotation.gtf> \
    -f <human_genome.fasta> \
    -gc <gnomad_constraint> \
    -out <output file prefix>
    -utr
```
## Output formats specification
### tab-separated (tsv) file
Each row represents annotation of a single variant in particular uORF (per uORF annotation).
#### Fields in the TSV output
1) #CHROM - contig name  
2) POS - position  
3) REF - reference allele
4) ALT - alternative allele
5) orf_start - uORF start position
6) orf_end - uORF end position
7) strand - strand dirrection defined as + (forward) or - (reverse).
8) gene_name - symnol gene ID from GTF annotation
9) transcript - transcript ID from uORF file
10) codon_type - ATG or non-ATG start codon in uORF 
11) overlapping_type - type of uORF (overlapping/non-overlapping) from input `.bed` file  
12) dist_from_orf_to_snp - distance from start of uORF to variant position  
13) consequence - type of effect 
14) main_cds_effect - effect of uORF variant on the main CDS (annotated for `frameshift` and `stop lost` variants only)
15) in_known_ORF - a binary flag indicating if a variant affects main ORF of the gene
16) pLI score (if gnomAD constraint annotation is provided)
17) LOEUF score (if gnomAD constraint annotation is provided)
18) INFO - old INFO field from inputed `.vcf` file  

### The Variant Call Format (vcf)
Add uBERT field to INFO fields of input vcf file.
#### FORMAT:
ORF_START|ORF_END|ORF_SYMB|ORF_CONSEQ|overlapped_CDS|utid|overlapping_type|dominance_type|codon_type
