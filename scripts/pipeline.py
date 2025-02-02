# pipeline.py

import logging
import pandas as pd
from typing import Optional, Dict, Tuple
from bedtools_utils import BedToolsWrapper
from coordinate_converter import CoordinateConverter
from variant_annotator import VariantAnnotator, Variant, ExonBoundary

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Pipeline:
    """Main pipeline for processing variants and converting coordinates."""
    
    def __init__(self, bed_file: str, vcf_file: str, gtf_file: str, fasta_file: str):
        """
        Initialize pipeline with input files.
        
        Args:
            bed_file: Path to BED file with uORF coordinates
            vcf_file: Path to VCF file with variants
            gtf_file: Path to GTF file with genomic features
            fasta_file: Path to FASTA file for sequence extraction
        """
        self.bed_file = bed_file
        self.vcf_file = vcf_file
        self.gtf_file = gtf_file
        self.fasta_file = fasta_file
        
        logger.info("Validating input files...")
        self._validate_input_files()
        
        self.bedtools = BedToolsWrapper()
        logger.info("Initializing CoordinateConverter...")
        self.converter = CoordinateConverter(bed_file, gtf_file)
        
        logger.info("Extracting sequences...")
        self.sequences = self.bedtools.get_fasta(bed_file, fasta_file)
        
        # Load BED file once
        self.bed_df = pd.read_csv(bed_file, sep='\t', header=None)

    def _validate_input_files(self):
        """Validate input file formats."""
        if not BedToolsWrapper.check_file_format(self.bed_file, 'bed'):
            raise ValueError(f"Invalid BED file: {self.bed_file}")
        if not BedToolsWrapper.check_file_format(self.vcf_file, 'vcf'):
            raise ValueError(f"Invalid VCF file: {self.vcf_file}")
        if not BedToolsWrapper.check_file_format(self.fasta_file, 'fasta'):
            raise ValueError(f"Invalid FASTA file: {self.fasta_file}")

    def _get_uorf_coordinates(self, transcript_id: str) -> Optional[Tuple[int, str]]:
        """Get uORF coordinates and strand from BED file."""
        try:
            uorf_row = self.bed_df[self.bed_df[3].str.contains(transcript_id)].iloc[0]
            return uorf_row[1], uorf_row[2], uorf_row[5]  # start, end, strand
        except (IndexError, KeyError):
            logger.warning(f"uORF coordinates not found for transcript {transcript_id}")
            return None

    def _get_transcript_info(self, transcript_id: str) -> Optional[dict]:
        """Get transcript information including sequences and coordinates."""
        # Get sequence from FASTA
        if transcript_id not in self.sequences:
            logger.warning(f"Sequence not found for transcript {transcript_id}")
            return None
        
        sequence = self.sequences[transcript_id]
        
        # Get transcript information
        transcript = self.converter.transcripts.get(transcript_id)
        if not transcript:
            logger.warning(f"Transcript information not found for {transcript_id}")
            return None

        # Get CDS coordinates
        cds_coords = self.converter.get_cds_coordinates(transcript_id)
        if not cds_coords:
            logger.warning(f"CDS coordinates not found for {transcript_id}")
            return None

        # Get uORF coordinates from BED file
        uorf_coords = self._get_uorf_coordinates(transcript_id)
        if not uorf_coords:
            return None
        
        uorf_start_genome, uorf_stop_genome, strand = uorf_coords
        
        # Convert genomic uORF coordinates to transcript coordinates
        uorf_start = self.converter.genome_to_transcript_pos(transcript_id, uorf_start_genome)
        uorf_stop = self.converter.genome_to_transcript_pos(transcript_id, uorf_stop_genome)
        
        if isinstance(uorf_start, str) or isinstance(uorf_stop, str):
            logger.warning(f"Failed to convert uORF coordinates for {transcript_id}")
            return None
        
        # Create ExonBoundary objects
        exon_boundaries = [
            ExonBoundary(
                start=exon.start,
                end=exon.start + exon.length - 1
            ) for exon in transcript.exons
        ]
        
        logger.info(f"Transcript info for {transcript_id}:")
        logger.info(f"  Sequence length: {len(sequence)}")
        logger.info(f"  uORF coords (genome): {uorf_start_genome}-{uorf_stop_genome}")
        logger.info(f"  uORF coords (transcript): {uorf_start}-{uorf_stop}")
        logger.info(f"  CDS coords (transcript): {cds_coords[0]}-{cds_coords[1]}")
        logger.info(f"  Strand: {strand}")
        
        return {
            'sequence': sequence,
            'exon_boundaries': exon_boundaries,
            'transcript': transcript,
            'uorf_start': uorf_start,
            'uorf_stop': uorf_stop,
            'cds_start': cds_coords[0],
            'cds_stop': cds_coords[1],
            'strand': strand
        }

    def process_variants(self) -> pd.DataFrame:
        """Process VCF variants and convert coordinates."""
        logger.info("Intersecting VCF with BED...")
        intersected = self.bedtools.intersect(self.vcf_file, self.bed_file)

        logger.info("Converting coordinates...")
        results = []
        total_rows = len(intersected)
        failed_conversions = 0

        for i, row in intersected.iterrows():
            if (i + 1) % 1000 == 0 or (i + 1) == total_rows:
                logger.info(f"Processed {i + 1}/{total_rows} variants")

            vcf_pos = row['pos']
            transcript_id = row['uorf_info'].split('|')[0].split('.')[0]

            logger.info(f"\nProcessing variant:")
            logger.info(f"  Transcript: {transcript_id}")
            logger.info(f"  Genomic position: {vcf_pos}")
            logger.info(f"  Reference: {row['ref']}")
            logger.info(f"  Alternative: {row['alt']}")

            if transcript_id not in self.converter.transcripts:
                logger.warning(f"Transcript {transcript_id} not found")
                continue

            # Convert genomic position to transcript position
            transcript_pos = self.converter.genome_to_transcript_pos(transcript_id, vcf_pos)
            logger.info(f"  Transcript position: {transcript_pos}")
            
            if transcript_pos == "NA":
                failed_conversions += 1
                continue

            # Get transcript information and sequence
            transcript_info = self._get_transcript_info(transcript_id)
            if not transcript_info:
                continue

            # Create Variant object for annotation
            variant = Variant(
                ref=row['ref'],
                alt=row['alt'],
                pos=vcf_pos,
                strand=transcript_info['strand']  # Используем strand из transcript_info
            )

            try:
                # Initialize annotator for this variant
                annotator = VariantAnnotator(
                    transcript_sequence=transcript_info['sequence'],
                    uorf_start=transcript_info['uorf_start'],
                    uorf_stop=transcript_info['uorf_stop'],
                    cds_start=transcript_info['cds_start'],
                    cds_stop=transcript_info['cds_stop'],
                    exon_boundaries=transcript_info['exon_boundaries']
                )
                
                # Annotate variant
                annotation_type, updated_exons = annotator.annotate_variant(
                    variant=variant,
                    transcript_pos=transcript_pos if isinstance(transcript_pos, int) else 0
                )
                annotation = annotation_type.value if annotation_type else 'Unknown'
                variant_type = variant.get_type().value
                
                logger.info(f"  Annotation result: {annotation}")
                logger.info(f"  Variant type: {variant_type}")
                
            except Exception as e:
                logger.warning(f"Failed to annotate variant at position {vcf_pos}: {str(e)}")
                annotation = 'Failed_to_annotate'
                variant_type = 'Unknown'

            # Convert transcript position back to genomic coordinates for validation
            recovered_genome_pos = self.converter.transcript_to_genome_pos(
                transcript_id, 
                transcript_pos
            )

            # Get nearest exon boundaries
            prev_exon_end, next_exon_start = self.converter.get_nearest_exon_boundaries(
                transcript_id, 
                vcf_pos
            )

            results.append({
                'Chromosome': row['chrom'],
                'Original_Genome_Position': vcf_pos,
                'Ref_Allele': row['ref'],
                'Alt_Allele': row['alt'],
                'Strand': transcript_info['strand'],
                'Transcript_ID': transcript_id,
                'Transcript_Position': transcript_pos,
                'Recovered_Genome_Position': recovered_genome_pos,
                'Coordinates_Match': str(vcf_pos) == str(recovered_genome_pos),
                'In_Exon': self.converter.check_position_in_exons(transcript_id, vcf_pos),
                'Nearest_Exon_Start': next_exon_start,
                'Nearest_Exon_End': prev_exon_end,
                'Variant_Type': variant_type,
                'Annotation': annotation
            })

        if failed_conversions > 0:
            logger.warning(
                f"Failed conversions: {failed_conversions}/{total_rows} "
                f"({(failed_conversions/total_rows)*100:.2f}%)"
            )

        return pd.DataFrame(results)

    def run(self, output_file: str) -> None:
        """
        Run the complete pipeline.
        
        Args:
            output_file: Path to output file
        """
        try:
            logger.info("Starting pipeline...")
            
            # Process variants
            results_df = self.process_variants()
            
            # Save results
            logger.info(f"Saving results to {output_file}")
            results_df.to_csv(output_file, sep='\t', index=False)
            
            # Log summary statistics
            if len(results_df) == 0:
                logger.warning("No variants were processed successfully")
                return
                
            total_variants = len(results_df)
            matching_coords = results_df['Coordinates_Match'].sum()
            in_exon_variants = results_df['In_Exon'].sum()
            
            logger.info("Pipeline Summary:")
            logger.info(f"Total variants processed: {total_variants}")
            logger.info(f"Variants with matching coordinates: {matching_coords} "
                       f"({(matching_coords/total_variants)*100:.2f}%)")
            logger.info(f"Variants in exons: {in_exon_variants} "
                       f"({(in_exon_variants/total_variants)*100:.2f}%)")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
            raise
        finally:
            logger.info("Pipeline finished.")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Process variants and convert coordinates')
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--fasta', required=True, help='Path to FASTA file')
    parser.add_argument('--output', required=True, help='Path to output file')

    args = parser.parse_args()

    try:
        pipeline = Pipeline(
            bed_file=args.bed,
            vcf_file=args.vcf,
            gtf_file=args.gtf,
            fasta_file=args.fasta
        )
        
        pipeline.run(output_file=args.output)
        
    except Exception as e:
        logger.error(f"Error running pipeline: {str(e)}", exc_info=True)
        raise