# uORF Variant Analysis Pipeline

A bioinformatics pipeline for analyzing the impact of genetic variants on upstream Open Reading Frames (uORFs) and predicting their consequences on the main Coding Sequence (CDS).

## Overview

This pipeline identifies and annotates variants that affect uORFs, which are small open reading frames located in the 5' untranslated regions (5' UTRs) of mRNAs. Variants in uORFs can impact translation regulation of the main protein-coding sequence, potentially leading to biological effects. The pipeline predicts consequences of variants on both the uORF itself and the downstream main CDS.

## Features

- Parses GTF and BED files to extract transcript and uORF information
- Analyzes variants in VCF format and intersects them with uORF regions
- Categorizes uORF variants into different consequence types (e.g., start loss, stop loss, frameshift)
- Predicts impact on main CDS (e.g., N-terminal extension, overlap truncation)
- Handles complex cases like boundary-crossing indels and splice variants
- Support for both canonical (ATG) and non-canonical start codons
- Generates comprehensive TSV output with annotation details
- Creates BED files for visualization in genome browsers

## Installation

### Prerequisites

- Python 3.7+
- Conda or Miniconda

### Setup Environment

1. Clone the repository:
```bash
git clone https://github.com/yourusername/uorf-variant-analysis.git
cd uorf-variant-analysis
```

2. Create and activate conda environment:
```bash
conda env create -f environment.yml
conda activate uorf-variant-analysis
```

## Usage

### Basic Usage

```bash
python pipeline.py --bed path/to/uorfs.bed --vcf path/to/variants.vcf \
                  --gtf path/to/annotation.gtf --fasta path/to/genome.fa \
                  --output-prefix results/output
```

### Command-line Options

| Option | Description |
|--------|-------------|
| `--bed` | Path to BED file with uORF coordinates |
| `--vcf` | Path to VCF file with variants |
| `--gtf` | Path to GTF annotation file |
| `--fasta` | Path to reference genome FASTA file |
| `--output-prefix` | Prefix for output files (.tsv and .bed will be appended) |
| `--uorf-type` | Filter uORFs by start codon type (`ALL`, `ATG`, or `NON-ATG`). Default: `ALL` |
| `--exclude-maincds-variants` | Exclude variants that are located within the main CDS region |
| `--debug` | Enable detailed debugging logs |

### Input File Requirements

#### BED File Format
The BED file should contain uORF coordinates with the following columns:
1. Chromosome
2. Start position (0-based)
3. End position (exclusive)
4. Name field (format: transcript_id|additional_info|start_codon_type)
5. Score (not used)
6. Strand

Example:
```
chr1    12345   12400   ENST00000123456|uORF1|ATG    0    +
```

#### VCF File
Standard VCF format with variants to be analyzed.

#### GTF File
Standard GTF format with gene annotations. The file must include CDS and exon features.

### Output Files

1. **TSV output** (output_prefix.tsv): Contains detailed annotation of each variant including:
   - Variant information (chromosome, position, rsID, alleles)
   - uORF information (coordinates, transcript ID)
   - Consequence on uORF (e.g., frameshift, start_lost)
   - Impact on main CDS (e.g., n_terminal_extension)
   - Codon changes

2. **BED output** (output_prefix.bed): Contains visualizable genomic regions affected by variants, which can be loaded into genome browsers.

## Core Components

### Module Overview

- **pipeline.py**: Main entry point for the pipeline
- **parsers.py**: Parsers for GTF and other genomic file formats
- **converters.py**: Handles conversion between genomic and transcript coordinates
- **processors.py**: Processes variants to determine their effects
- **annotator.py**: Annotates variants with biological consequences
- **models.py**: Data structures for genomic and transcript features
- **transcript_sequence.py**: Handles transcript sequences and uORF extraction

### Key Classes

- **Pipeline**: Main controller that orchestrates the analysis workflow
- **CoordinateConverter**: Converts between genomic and transcript coordinates
- **VariantProcessor**: Processes variants and determines their effects
- **VariantAnnotator**: Annotates variants with biological consequences
- **Transcript**: Represents transcript data with coordinate mappings
- **TranscriptSequence**: Handles sequence extraction and manipulation

## Variant Classification

### uORF Consequences

The pipeline classifies variants into the following consequence types:

- **START_LOST**: Loss of uORF start codon
- **STOP_LOST**: Loss of uORF stop codon
- **STOP_GAINED**: Creation of a premature stop codon
- **FRAMESHIFT**: Indel causing a shift in reading frame
- **DELETION_AND_STOP_LOST**: Complex cases where both deletion and stop loss occur
- **MISSENSE**: Nonsynonymous variants changing amino acid
- **SYNONYMOUS**: Synonymous variants preserving amino acid
- **SPLICE_SITE**: Variants affecting splice sites
- **INFRAME_DELETION**: Deletions that maintain the reading frame
- **INFRAME_INSERTION**: Insertions that maintain the reading frame

### Main CDS Impact

Predicts how the uORF variant affects the main CDS:

- **N_TERMINAL_EXTENSION**: Extension of protein N-terminus
- **OUT_OF_FRAME_OVERLAP**: Out-of-frame overlap with main CDS
- **UORF_PRODUCT_TRUNCATION**: Truncation of uORF product
- **UORF_PRODUCT_EXTENSION**: Extension of uORF product
- **STOP_GAINED**: Introduction of premature stop codon
- **OVERLAP_EXTENSION**: Extension of uORF-CDS overlap
- **OVERLAP_TRUNCATION**: Truncation of uORF-CDS overlap
- **OVERLAP_ELIMINATION**: Elimination of uORF-CDS overlap
- **MAIN_CDS_UNAFFECTED**: No effect on main CDS

## Examples

### Example 1: Analyze variants with default settings

```bash
python pipeline.py --bed data/uorfs.bed --vcf data/variants.vcf \
                  --gtf data/gencode.v38.annotation.gtf --fasta data/GRCh38.p13.genome.fa \
                  --output-prefix results/all_variants
```

### Example 2: Analyze only variants in ATG-start uORFs

```bash
python pipeline.py --bed data/uorfs.bed --vcf data/variants.vcf \
                  --gtf data/gencode.v38.annotation.gtf --fasta data/GRCh38.p13.genome.fa \
                  --output-prefix results/atg_uorf_variants --uorf-type ATG
```

### Example 3: Exclude variants that fall within main CDS

```bash
python pipeline.py --bed data/uorfs.bed --vcf data/variants.vcf \
                  --gtf data/gencode.v38.annotation.gtf --fasta data/GRCh38.p13.genome.fa \
                  --output-prefix results/non_cds_variants --exclude-maincds-variants
```

## Cite

If you use this pipeline in your research, please cite:

[PMID: 36651276](https://pubmed.ncbi.nlm.nih.gov/36651276/)
