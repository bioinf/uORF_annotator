import logging
import argparse
import pandas as pd
from pybedtools import BedTool
import pysam

from converters import CoordinateConverter
from processors import VariantProcessor


class Pipeline:
    """Main pipeline for variant analysis and annotation."""
    
    def __init__(self, bed_file: str, vcf_file: str, gtf_file: str, fasta_file: str):
        """Initialize pipeline with input files."""
        self.bed_file = bed_file
        self.vcf_file = vcf_file
        self.fasta_file = fasta_file
        self.fasta = pysam.FastaFile(fasta_file)
        self.converter = CoordinateConverter(bed_file, gtf_file)
        self.processor = VariantProcessor(self.converter, self.fasta)

    def process_variants(self) -> pd.DataFrame:
        """Process and annotate variants."""
        intersected = self._intersect_files()
        if intersected.empty:
            logging.error("No intersections found between VCF and BED!")
            return pd.DataFrame()

        results = []
        total_rows = len(intersected)

        for i in range(total_rows):
            if (i >= 0 and i % 1000 == 0) or (i == total_rows - 1):
                logging.info(f"Processed {i}/{total_rows} variants")

            # Process variant now returns a list of results for all variant-uORF pairs
            variant_results = self.processor.process_variant(intersected.iloc[i])
            if variant_results:
                # Add all results instead of just the first one
                results.extend(variant_results)

        results_df = pd.DataFrame(results)
        
        # Sort results by chromosome and position (removed uORF_ID from sorting)
        if not results_df.empty:
            results_df = results_df.sort_values(['Chromosome', 'Original_Genome_Position'])
        
        return results_df

    def _intersect_files(self) -> pd.DataFrame:
        """Intersect VCF and BED files using bedtools."""
        vcf = BedTool(self.vcf_file)
        bed = BedTool(self.bed_file)
        intersection = vcf.intersect(bed, wa=True, wb=True)
        
        df = pd.read_csv(intersection.fn, sep='\t', header=None)
        df.columns = [f'col{i}' for i in range(len(df.columns))]
        return df


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description='Analyze variants in uORFs and predict their impact on main CDS'
    )
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('--output', required=True, help='Path to output file')

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    try:
        logging.info("Initializing pipeline...")
        pipeline = Pipeline(args.bed, args.vcf, args.gtf, args.fasta)
        
        logging.info("Processing variants...")
        results = pipeline.process_variants()
        
        if not results.empty:
            results.to_csv(args.output, sep='\t', index=False)
            logging.info(f"Results successfully saved to {args.output}")
            logging.info(f"Total variant-uORF pairs processed: {len(results)}")
        else:
            logging.warning("No variants were processed successfully.")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}", exc_info=True)
        raise


if __name__ == '__main__':
    main()