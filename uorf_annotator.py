import logging
import argparse
import pandas as pd
from pybedtools import BedTool
import pysam
import os

from scripts.converters import CoordinateConverter
from scripts.processors import VariantProcessor
from scripts.annotator import VariantAnnotator


class Pipeline:
    """Main pipeline for variant analysis and annotation."""
    
    def __init__(self, bed_file: str, vcf_file: str, gtf_file: str, fasta_file: str, 
                output_prefix: str, uorf_type: str = "ALL", exclude_maincds_variants: bool = False,
                debug_mode: bool = False):
        """Initialize pipeline with input files."""
        self.bed_file = bed_file
        self.vcf_file = vcf_file
        self.fasta_file = fasta_file
        self.output_prefix = output_prefix
        self.debug_mode = debug_mode
        self.uorf_type = uorf_type.upper()  # Store the uORF type option
        self.exclude_maincds_variants = exclude_maincds_variants  # Store the flag for exclusion
        
        # Validate uORF type value
        if self.uorf_type not in ["ALL", "ATG", "NON-ATG"]:
            logging.warning(f"Invalid uORF type '{self.uorf_type}', defaulting to 'ALL'")
            self.uorf_type = "ALL"
        
        # Generate output filenames - ensure we don't double-add extensions
        self.tsv_output = f"{output_prefix}.tsv"
        self.bed_output = f"{output_prefix}.bed"
        
        # Initialize fasta reference and converter before creating output files
        self.fasta = pysam.FastaFile(fasta_file)
        
        # Use a copy of the bed file for analysis to avoid file handle issues
        self.converter = CoordinateConverter(bed_file, gtf_file, uorf_type=self.uorf_type, debug_mode=debug_mode)
        
        # Initialize the processor with the output BED file path and the exclusion flag
        self.processor = VariantProcessor(self.converter, self.fasta, bed_output=self.bed_output, 
                                        exclude_maincds_variants=self.exclude_maincds_variants,
                                        debug_mode=debug_mode)

    def process_variants(self) -> pd.DataFrame:
        """Process and annotate variants."""
        # Create an empty BED output file right before processing
        # This ensures any file handles from reading are closed
        with open(self.bed_output, 'w') as bed_file:
            # Just create an empty file - will be appended to during processing
            pass
            
        try:
            intersected = self._intersect_files()
            if intersected.empty:
                logging.warning("No intersections found between VCF and BED!")
                return pd.DataFrame()
        except Exception as e:
            logging.error(f"Error during file intersection: {str(e)}")
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
        logging.info(f"Using uORF type filter: {self.uorf_type}")

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
            results_df = results_df.sort_values(['Chromosome', 'Variant_Genomic_Position'])
        
        return results_df

    def _intersect_files(self) -> pd.DataFrame:
        """Intersect VCF and BED files using bedtools."""
        try:
            # Add the -wb option to include all columns from the VCF file including INFO field
            vcf = BedTool(self.vcf_file)
            bed = BedTool(self.bed_file)
            
            # Change intersection to include all VCF fields
            intersection = vcf.intersect(bed, wb=True)
            
            # Check if the intersection produced any results
            try:
                # Try to read the first line to see if there are any results
                with open(intersection.fn, 'r') as f:
                    first_line = f.readline().strip()
                if not first_line:
                    logging.warning("Intersection result is empty")
                    return pd.DataFrame()  # Return empty DataFrame
            except:
                logging.warning("Failed to read intersection results, returning empty DataFrame")
                return pd.DataFrame()  # Return empty DataFrame if can't read
            
            # Try to read the intersection results into a DataFrame
            try:
                df = pd.read_csv(intersection.fn, sep='\t', header=None)
                df.columns = [f'col{i}' for i in range(len(df.columns))]
                
                if self.debug_mode:
                    logging.debug(f"Intersection result: first few rows:\n{df.head()}")
                    logging.debug(f"Columns: {df.columns}")
                
                return df
            except pd.errors.EmptyDataError:
                logging.warning("Intersection result is empty (pandas EmptyDataError)")
                return pd.DataFrame()  # Return empty DataFrame
            
        except Exception as e:
            logging.error(f"Error during intersection: {str(e)}")
            return pd.DataFrame()  # Return empty DataFrame on any error

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
            # Create empty output files so the pipeline completes cleanly
            with open(self.tsv_output, 'w') as f:
                f.write("No variants were processed successfully\n")
            # BED file was already created at the beginning of process_variants()
            logging.info(f"Empty result files created at {self.tsv_output} and {self.bed_output}")


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
    
    # Add uORF type filter argument
    parser.add_argument('--uorf-type', default='ALL', choices=['ALL', 'ATG', 'NON-ATG'], 
                       help='Filter uORFs by start codon type (ALL, ATG, or NON-ATG). Default: ALL')
    parser.add_argument('--exclude-maincds-variants', action='store_true', 
                       help='Exclude variants that are located within the main CDS region')
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
                        uorf_type=args.uorf_type, exclude_maincds_variants=args.exclude_maincds_variants,
                        debug_mode=args.debug)
        
        logging.info("Processing variants...")
        results = pipeline.process_variants()

        # Save results to both TSV and BED files
        pipeline.save_results(results)
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}", exc_info=True)
        raise


if __name__ == '__main__':
    main()