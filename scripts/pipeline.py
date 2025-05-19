import logging
import argparse
import pandas as pd
from pybedtools import BedTool
import pysam
import os

from converters import CoordinateConverter
from processors import VariantProcessor
from annotator import VariantAnnotator


class Pipeline:
    """Main pipeline for variant analysis and annotation."""
    
    def __init__(self, bed_file: str, vcf_file: str, gtf_file: str, fasta_file: str, 
                 output_prefix: str, frame_filter: str = "ALL", debug_mode: bool = False):
        """Initialize pipeline with input files."""
        self.bed_file = bed_file
        self.vcf_file = vcf_file
        self.fasta_file = fasta_file
        self.output_prefix = output_prefix
        self.debug_mode = debug_mode
        self.frame_filter = frame_filter.upper()  # Store the frame filter option
        
        # Validate frame filter value
        if self.frame_filter not in ["ALL", "ATG", "NON-ATG"]:
            logging.warning(f"Invalid frame filter '{self.frame_filter}', defaulting to 'ALL'")
            self.frame_filter = "ALL"
        
        # Generate output filenames - ensure we don't double-add extensions
        self.tsv_output = f"{output_prefix}.tsv"
        self.bed_output = f"{output_prefix}.bed"
        
        # Initialize fasta reference and converter before creating output files
        self.fasta = pysam.FastaFile(fasta_file)
        
        # Use a copy of the bed file for analysis to avoid file handle issues
        self.converter = CoordinateConverter(bed_file, gtf_file, frame_filter=self.frame_filter, debug_mode=debug_mode)
        
        # Initialize the processor with the output BED file path
        # but don't create the file until we're ready to write to it
        self.processor = VariantProcessor(self.converter, self.fasta, bed_output=self.bed_output, debug_mode=debug_mode)

    def process_variants(self) -> pd.DataFrame:
        """Process and annotate variants."""
        # Create an empty BED output file right before processing
        # This ensures any file handles from reading are closed
        with open(self.bed_output, 'w') as bed_file:
            # Just create an empty file - will be appended to during processing
            pass
            
        intersected = self._intersect_files()
        if intersected.empty:
            logging.error("No intersections found between VCF and BED!")
            return pd.DataFrame()

        results = []
        # Keep track of which variants have been processed to avoid duplicates
        processed_variants = set()
        
        # Group the intersection results by VCF positions
        # This allows us to process each variant only once, considering all its BED overlaps
        position_groups = {}
        for idx, row in intersected.iterrows():
            variant_key = f"{row['col0']}:{row['col1']}:{row['col3']}>{row['col4']}"
            if variant_key not in position_groups:
                position_groups[variant_key] = []
            position_groups[variant_key].append(idx)
        
        logging.info(f"Found {len(position_groups)} unique variants after grouping")
        logging.info(f"Using frame filter: {self.frame_filter}")

        # Process each unique variant position
        for variant_key, row_indices in position_groups.items():
            chrom, pos, ref_alt = variant_key.split(':', 2)
            logging.info(f"Processing variant {chrom}:{pos} {ref_alt}")
            
            # Collect all BED intersections for this variant
            combined_row = None
            bed_entries = []
            
            for idx in row_indices:
                row = intersected.iloc[idx]
                if combined_row is None:
                    # Use the first row as the base
                    combined_row = row.copy()
                
                # Store BED information
                bed_entry = {
                    'start': int(row['col9']),
                    'end': int(row['col10']),
                    'name': row['col11'],
                    'strand': row['col13']
                }
                bed_entries.append(bed_entry)
            
            if combined_row is None:
                logging.warning(f"No valid data found for variant {variant_key}")
                continue
                
            # Create a marker for the processor to use all BED entries
            combined_row['all_bed_entries'] = bed_entries
            
            # Process the combined variant information
            variant_results = self.processor.process_variant(combined_row)
            if variant_results:
                # Add the results and mark this variant as processed
                results.extend(variant_results)
                processed_variants.add(variant_key)
                logging.info(f"Added {len(variant_results)} results for {variant_key}")
            else:
                logging.warning(f"No results generated for {variant_key}")

        logging.info(f"Processed {len(processed_variants)} unique variants, generated {len(results)} total results")
        
        results_df = pd.DataFrame(results)
        
        # Sort results by chromosome and position (removed uORF_ID from sorting)
        if not results_df.empty:
            results_df = results_df.sort_values(['Chromosome', 'Original_Genome_Position'])
        
        return results_df

    # In the _intersect_files method:
    def _intersect_files(self) -> pd.DataFrame:
        """Intersect VCF and BED files using bedtools."""
        # Add the -wb option to include all columns from the VCF file including INFO field
        vcf = BedTool(self.vcf_file)
        bed = BedTool(self.bed_file)
        # Change intersection to include all VCF fields
        intersection = vcf.intersect(bed, wb=True)
        
        df = pd.read_csv(intersection.fn, sep='\t', header=None)
        df.columns = [f'col{i}' for i in range(len(df.columns))]
        
        if self.debug_mode:
            logging.debug(f"Intersection result: first few rows:\n{df.head()}")
            logging.debug(f"Columns: {df.columns}")
        
        return df

    def save_results(self, results_df: pd.DataFrame) -> None:
        """Save results to output files."""
        if not results_df.empty:
            # Save TSV results
            results_df.to_csv(self.tsv_output, sep='\t', index=False)
            logging.info(f"Results successfully saved to {self.tsv_output}")
            logging.info(f"BED file generated at {self.bed_output}")
            logging.info(f"Total variant-uORF pairs processed: {len(results_df)}")
        else:
            logging.warning("No variants were processed successfully.")


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description='Analyze variants in uORFs and predict their impact on main CDS'
    )
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--fasta', required=True, help='Path to reference FASTA')
    
    # Support both old and new output arguments for backward compatibility
    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument('--output', help='Path to output file (legacy mode)')
    output_group.add_argument('--output-prefix', help='Prefix for output files (.tsv and .bed will be appended)')
    
    # Add frame filter argument
    parser.add_argument('--frame-filter', default='ALL', choices=['ALL', 'ATG', 'NON-ATG'], 
                       help='Filter uORFs by start codon type (ALL, ATG, or NON-ATG). Default: ALL')
    
    # Add debug mode argument
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')

    args = parser.parse_args()

    # Set logging level based on the debug flag
    log_level = logging.DEBUG if args.debug else logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    try:
        # Determine output prefix from either argument
        output_prefix = args.output_prefix if args.output_prefix else args.output
        
        # Clean up the output prefix to remove any extensions
        # This ensures we don't end up with duplicate extensions
        if output_prefix.endswith('.tsv'):
            output_prefix = output_prefix[:-4]
        elif output_prefix.endswith('.bed'):
            output_prefix = output_prefix[:-4]
            
        logging.info("Initializing pipeline...")
        pipeline = Pipeline(args.bed, args.vcf, args.gtf, args.fasta, output_prefix, 
                           frame_filter=args.frame_filter, debug_mode=args.debug)
        
        logging.info("Processing variants...")
        results = pipeline.process_variants()

        # Save results to both TSV and BED files
        pipeline.save_results(results)
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}", exc_info=True)
        raise


if __name__ == '__main__':
    main()