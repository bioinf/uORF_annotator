# uORF_annotator
uORF_annotator is the tool for annotation of upstream translation initiation sites in upstream open reading frames (uORFs) of human genome, which predicted by [uBert model](https://github.com/skoblov-lab/uBERTa).
## Conda environment
Install all dependencies from `requirements.yml` as new conda environment.
```
conda env create -f requirements.yml
```
## Required input data
* VCF file of variants for further annotation
* BED file with available uORFs
* GTF file with genomic features annotation (for example [GENCODE](https://www.gencodegenes.org/human/))
## Run example
```
python uORF_annotator.py \
    -i <input_variants.vcf> \
    -b <input_uORFs.bed> \
    -g <input_annotation.gtf> \
    -f <human_genome.fasta> \
    -ot <output.tsv> \
    -ov <output.vsc>
```
## Output formats specification
### tab-separated (tsv) file
Each row represents one variant in particular uORF (per uORF annotation)
#### Fields
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
14) overlapping_type - type of overlapping from input `.bed` file  
15) codon_type - ATG or non-ATG type  
16) dist_from_orf_to_snp - distance from begin of uORF to variant position  
17) symbol - schematic indication of the variant event  
18) consequence - type of effect 
19) in_known_ORF - is variant also in main ORF
20) in_known_CDS - is variant also in main CDS
### The Variant Call Format (vcf)
Add uBERT field to INFO fields of input vcf file.
#### FORMAT:
ORF_START|ORF_END|ORF_SYMB|ORF_CONSEQ|overlapped_CDS|utid|overlapping_type|codon_type
