# uORF Annotator v. 1.0
*uORF Annotator* is the tool for annotation of the functional impact of genetic variants in upstream open reading frames (uORFs) in the human genome, which were manually annotated based on publicly available Ribo-seq and other data types in 3641 OMIM genes.

## Conda environment
Install all dependencies from `requirements.yml` as new conda environment.
```
conda env create -f requirements.yml
```
## Required input data
* VCF file of variants for further annotation
* BED file with available uORFs (`sorted.v4.bed` in this repository)
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

Two TSV outputs are generated - one for ATG-started uORFs and one - for non-ATG-started ones. Each row represents annotation of a single variant in particular uORF (per uORF annotation). Fields in the file have the following content:

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

### Varinat call format (VCF)

The generated VCF output contains all variants affecting uORF sequences. Each variant is annotated with the following INFO fields: `uORFs`, `uORFs_ATG`, `uORFs_eff`. The description of fields is given below:

* `uORFs` - a full consequence annotation for each variant-uORF combination. Format: 'ORF_START|ORF_END|ORF_SYMB|ORF_CONSEQ|main_cds_effect|in_known_CDS|in_known_ORF|utid|overlapping_type|dominance_type|codon_type'
* `uORFs_ATG` - a flag indicating if a variant falls within at least one ATG-starting uORF.
* `uORFs_eff` - a short notation of how a change in uORF structure resulting from a variant affects the main coding part (СDS) of a gene. ext - N-terminal extension, overl - out-of-frame overlap, activ - overlap removal with possible main ORF activation, unaff - no effect on main CDS. If the variant falls into more than one uORF, the effects on them are listed through &

### BED format

*uORF Annotator* generates two BED files with uORFs affected by variants that alter the uORF length (one file contains ATG uORFs and the other contains non-ATG-started uORFs). Both BED files contain two entries for each affected uORF: 
1) initial uORF, its `name` field format: uORF_unique_number-gene_name|uORF_type|start_codon_type(ATG/non-AT), filled with black color; 
2) resulting uORF after introduction of a variant,  its `name` field format: uORF_unique_number-gene_name|variant|variant_type|main_CDS_effect,  filled with different colors depending on the effect. Color legend: 
* Grey features - cases when the variant does not change the overlap between uORF and main CDS.
* Orange features - cases when (a) uORF-truncating variant eliminates the existing overlap between uORF and main CDS; or (b) variant leads to the production of a chimeric protein product of the gene, possessing an extension at the N-terminus resulting from uORF translation
* Red features - cases where variant leads to the appearance of a new overlapping segment between uORF and main gene CDS, with the two sequences translated in different frames.


## Supplementary data files

This repository contains two additional files:
1) `Annotated_uORFs_and_alt.CDS.starts_v4.8.bed` - Manually annotated alternative open reading frames (including non-overlapping and overlapping uORFs, CDS_extensions and CDS_truncations) found in 3641 human genes from the OMIM database . BED-file, field "name" contains information about 'gene_name|isoform|type_of_ORF|start_codon'.
2) `High-confidence_uORFs_v2.bed` - List of high confidence uORFs found in 3641 human genes from the OMIM database. In this list, we included only uORFs predicted in at least two out of four different studies (the present study, Ji et al. (PMID: 26687005), McGillivray et al. (PMID: 29562350), Scholz et al. (PMID: 31513641)). BED-file, field "name" contains information about 'gene_name|isoform|type_of_uORF|type_of_start_codon(ATG/non-ATG)'.
