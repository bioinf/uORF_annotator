import logging
from typing import Optional, Dict
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact
from transcript_sequence import TranscriptSequence

class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile):
        """Initialize the variant processor with coordinate converter and reference sequence."""
        self.converter = converter
        self.fasta = fasta
        self.annotator = None
        self.debug_info = {}

    def process_variant(self, row: pd.Series) -> Optional[Dict]:
        """Process a variant to determine its effect on uORF and main CDS."""
        try:
            # Extract variant information
            chrom = row['col0']
            vcf_pos = int(row['col1'])
            ref_allele = row['col3']
            alt_allele = row['col4']
            
            logging.debug(f"Processing variant at {chrom}:{vcf_pos} {ref_allele}>{alt_allele}")
            
            # Get transcript information
            original_transcript_id = self._extract_transcript_id(row['col11'])
            
            # Find all matching transcripts, including those with multiple uORFs
            matching_transcripts = [
                tid for tid in self.converter.transcripts.keys() 
                if tid == original_transcript_id or tid.startswith(f"{original_transcript_id}_uorf_")
            ]
            
            if not matching_transcripts:
                logging.warning(f"No transcripts found for {original_transcript_id}")
                return None
            
            logging.debug(f"Found {len(matching_transcripts)} matching transcripts: {matching_transcripts}")
            
            # Process each matching transcript-uORF pair
            results = []
            
            for transcript_id in matching_transcripts:
                logging.debug(f"Processing transcript-uORF: {transcript_id}")
                
                transcript_obj = self.converter.transcripts[transcript_id]
                
                # Log debug information
                self._log_debug_info(transcript_id, transcript_obj, vcf_pos)

                # Check if transcript was extended to include uORF
                was_extended = getattr(transcript_obj, 'was_extended', False)
                
                # Get transcript position for the variant
                variant_coords = transcript_obj.get_coordinates(vcf_pos)
                if not variant_coords:
                    logging.debug(f"Position {vcf_pos} not in transcript coordinates for {transcript_id}, skipping")
                    continue  # Skip to next transcript
                    
                # Check if uORF coordinates are still missing after potential extension
                if transcript_obj.uorf_start is None or transcript_obj.uorf_end is None:
                    logging.debug(f"Missing uORF transcript coordinates for {transcript_id}, skipping")
                    continue  # Skip to next transcript
                    
                # Create TranscriptSequence object
                transcript_seq = TranscriptSequence(transcript_obj, self.fasta, chrom)
                if not transcript_seq.sequence:
                    logging.debug(f"Could not extract transcript sequence for {transcript_id}, skipping")
                    continue
                    
                # Check if uORF region was extracted successfully
                if not transcript_seq.uorf_region:
                    logging.debug(f"Failed to extract uORF region for {transcript_id}, skipping")
                    continue
                    
                # Initialize annotator with TranscriptSequence
                self.annotator = VariantAnnotator(transcript_seq)

                # Check for overlap between uORF and mainCDS
                overlaps_maincds = False
                if (transcript_obj.uorf_end is not None and 
                    transcript_obj.mainorf_start is not None):
                    overlaps_maincds = self.annotator.does_overlap_maincds(
                        transcript_obj.uorf_end, transcript_obj.mainorf_start
                    )

                # Get codon change
                codon_change = self.annotator.get_codon_change({
                    'position': variant_coords.transcript,
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele
                })
                    
                # Prepare variant data with all information
                variant_data = {
                    'position': variant_coords.transcript,
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele,
                    'strand': transcript_obj.strand,
                    'uorf_start': transcript_obj.uorf_start,
                    'uorf_end': transcript_obj.uorf_end,
                    'maincds_start': transcript_obj.mainorf_start,
                    'maincds_end': transcript_obj.mainorf_end,
                    'uorf_start_genomic': transcript_obj.uorf_start_genomic,
                    'uorf_end_genomic': transcript_obj.uorf_end_genomic,
                    'maincds_start_genomic': transcript_obj.mainorf_start_genomic,
                    'maincds_end_genomic': transcript_obj.mainorf_end_genomic,
                    'codon_change': codon_change,
                    'overlaps_maincds': overlaps_maincds,
                    'transcript_extended': was_extended
                }

                # Get consequences and impacts
                uorf_consequence = self.annotator.get_consequence(variant_data)
                maincds_impact = None
                if uorf_consequence:
                    maincds_impact = self.annotator.predict_impact(variant_data, uorf_consequence)
                
                # For output, we need to ensure biological orientation (5' to 3')
                # We use the original coordinates (before any internal fixing)
                # to maintain the proper biological meaning
                uorf_start_transcript = getattr(transcript_seq, 'original_uorf_start', transcript_obj.uorf_start)
                uorf_end_transcript = getattr(transcript_seq, 'original_uorf_end', transcript_obj.uorf_end)
                maincds_start_transcript = getattr(transcript_seq, 'original_mainorf_start', transcript_obj.mainorf_start)
                maincds_end_transcript = getattr(transcript_seq, 'original_mainorf_end', transcript_obj.mainorf_end)
                
                # For negative strand, ensure coordinates reflect biological orientation (5' to 3')
                if transcript_obj.strand == '-':
                    # Only swap if they weren't swapped before
                    if uorf_start_transcript < uorf_end_transcript:
                        uorf_start_transcript, uorf_end_transcript = uorf_end_transcript, uorf_start_transcript
                    
                    if maincds_start_transcript < maincds_end_transcript:
                        maincds_start_transcript, maincds_end_transcript = maincds_end_transcript, maincds_start_transcript
                    
                # Extract the original transcript ID if this is a derived uORF transcript
                display_transcript_id = original_transcript_id
                    
                result = {
                    'Chromosome': chrom,
                    'Original_Genome_Position': vcf_pos,
                    'Transcript_Position': variant_coords.transcript,
                    'Ref_Allele': ref_allele,
                    'Alt_Allele': alt_allele,
                    'Strand': transcript_obj.strand,
                    'Transcript_ID': display_transcript_id,
                    'uORF_Start_Genomic': transcript_obj.uorf_start_genomic,
                    'uORF_End_Genomic': transcript_obj.uorf_end_genomic,
                    'uORF_Start_Transcript': uorf_start_transcript,
                    'uORF_End_Transcript': uorf_end_transcript,
                    'mainCDS_Start_Genomic': transcript_obj.mainorf_start_genomic,
                    'mainCDS_End_Genomic': transcript_obj.mainorf_end_genomic,
                    'mainCDS_Start_Transcript': maincds_start_transcript,
                    'mainCDS_End_Transcript': maincds_end_transcript,
                    'Codon_Change': codon_change,
                    'uORF_Consequence': uorf_consequence.value if uorf_consequence else 'None',
                    'uORF_mainCDS_Overlap': 'overlapping' if overlaps_maincds else 'non_overlapping',
                    'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None'
                }
                
                # Add additional information if transcript was extended
                if was_extended:
                    result['Transcript_Extended'] = 'Yes'
                    
                results.append(result)
                logging.debug(f"Successfully processed variant for transcript-uORF: {transcript_id}")
            
            # Return the first successful result, or None if all failed
            if results:
                return results[0]  # For backward compatibility, return first result
            return None
                    
        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}")
            return None

    def _log_debug_info(self, transcript_id, transcript_obj, vcf_pos):
        """Log debug information about transcript coordinates."""
        debug_info = {
            'transcript_id': transcript_id,
            'strand': transcript_obj.strand,
            'uorf_start_g': transcript_obj.uorf_start_genomic,
            'uorf_end_g': transcript_obj.uorf_end_genomic,
            'uorf_start_t': transcript_obj.uorf_start,
            'uorf_end_t': transcript_obj.uorf_end,
            'variant_pos_g': vcf_pos,
            'was_extended': getattr(transcript_obj, 'was_extended', False)
        }
        
        # Store debug info
        self.debug_info[transcript_id] = debug_info
        
        # Check for potential issues
        if transcript_obj.uorf_start is None and transcript_obj.uorf_start_genomic is not None:
            logging.warning(f"uORF genomic coordinates exist but transcript coordinates are None for {transcript_id}")
            
            # Show transcript boundaries
            if transcript_obj.genome_to_transcript:
                min_pos = min(transcript_obj.genome_to_transcript.keys())
                max_pos = max(transcript_obj.genome_to_transcript.keys())
                
                # Check if uORF is outside transcript bounds
                if transcript_obj.uorf_start_genomic < min_pos or transcript_obj.uorf_start_genomic > max_pos:
                    logging.warning(f"uORF start ({transcript_obj.uorf_start_genomic}) is outside transcript bounds ({min_pos}-{max_pos})")
                if transcript_obj.uorf_end_genomic < min_pos or transcript_obj.uorf_end_genomic > max_pos:
                    logging.warning(f"uORF end ({transcript_obj.uorf_end_genomic}) is outside transcript bounds ({min_pos}-{max_pos})")
        
        if transcript_obj.uorf_start is not None and transcript_obj.uorf_end is not None:
            if transcript_obj.strand == '+' and transcript_obj.uorf_start > transcript_obj.uorf_end:
                logging.warning(f"Suspicious positive strand coordinates: uORF_start > uORF_end for {transcript_id}")
            elif transcript_obj.strand == '-' and transcript_obj.uorf_start < transcript_obj.uorf_end:
                logging.warning(f"Suspicious negative strand coordinates: uORF_start < uORF_end for {transcript_id}")
            
            # Check for invalid coordinates
            if transcript_obj.uorf_start < 1 or transcript_obj.uorf_end < 1:
                logging.warning(f"Invalid transcript coordinates: less than 1 for {transcript_id}")

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        """Extract transcript ID from BED name field."""
        return bed_name.split('|')[0].split('.')[0]