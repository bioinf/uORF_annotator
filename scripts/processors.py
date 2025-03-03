import logging
from typing import Optional, Dict
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact
from transcript_sequence import TranscriptSequence

class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile):
        """
        Initialize the variant processor with coordinate converter and reference sequence.
        
        Args:
            converter: CoordinateConverter object for genomic/transcript coordinate conversion
            fasta: FastaFile object containing reference sequences
        """
        self.converter = converter
        self.fasta = fasta
        self.annotator = None
        self.debug_info = {}

    def process_variant(self, row: pd.Series) -> Optional[Dict]:
        """
        Process a variant to determine its effect on uORF and main CDS.
        
        Args:
            row: Pandas Series containing variant information
            
        Returns:
            Dictionary with annotated variant information, or None if processing fails
        """
        try:
            # Extract variant information
            chrom = row['col0']
            vcf_pos = int(row['col1'])
            ref_allele = row['col3']
            alt_allele = row['col4']
            
            # Get transcript information
            transcript_id = self._extract_transcript_id(row['col11'])
            matching_transcripts = [
                tid for tid in self.converter.transcripts.keys() 
                if tid == transcript_id
            ]

            if not matching_transcripts:
                logging.warning(f"Transcript {transcript_id} not found")
                return None

            matched_transcript = matching_transcripts[0]
            transcript_obj = self.converter.transcripts[matched_transcript]

            # Log debug information
            self._log_debug_info(transcript_id, transcript_obj, vcf_pos)

            # Get transcript position
            variant_coords = transcript_obj.get_coordinates(vcf_pos)
            if not variant_coords:
                logging.warning(f"Could not convert position {vcf_pos} to transcript coordinates")
                return None

            if transcript_obj.uorf_start is None or transcript_obj.uorf_end is None:
                logging.warning(f"Missing uORF transcript coordinates for {transcript_id}")
                return None

            # Create TranscriptSequence object
            transcript_seq = TranscriptSequence(transcript_obj, self.fasta, chrom)
            if not transcript_seq.sequence or not transcript_seq.uorf_region:
                logging.warning(f"Could not extract transcript sequence for {transcript_id}")
                return None
                    
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
                'overlaps_maincds': overlaps_maincds
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
                
            return {
                'Chromosome': chrom,
                'Original_Genome_Position': vcf_pos,
                'Transcript_Position': variant_coords.transcript,
                'Ref_Allele': ref_allele,
                'Alt_Allele': alt_allele,
                'Strand': transcript_obj.strand,
                'Transcript_ID': matched_transcript,
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
                    
        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}")
            return None

    def _log_debug_info(self, transcript_id, transcript_obj, vcf_pos):
        """
        Log debug information about transcript coordinates.
        
        Args:
            transcript_id: ID of the transcript
            transcript_obj: Transcript object
            vcf_pos: Variant position in genome coordinates
        """
        debug_info = {
            'transcript_id': transcript_id,
            'strand': transcript_obj.strand,
            'uorf_start_g': transcript_obj.uorf_start_genomic,
            'uorf_end_g': transcript_obj.uorf_end_genomic,
            'uorf_start_t': transcript_obj.uorf_start,
            'uorf_end_t': transcript_obj.uorf_end,
            'variant_pos_g': vcf_pos
        }
        
        # Store debug info
        self.debug_info[transcript_id] = debug_info
        
        # Check for potential issues
        if transcript_obj.uorf_start is not None and transcript_obj.uorf_end is not None:
            if transcript_obj.strand == '+' and transcript_obj.uorf_start > transcript_obj.uorf_end:
                logging.warning(f"Suspicious positive strand coordinates: uORF_start > uORF_end for {transcript_id}")
            elif transcript_obj.strand == '-' and transcript_obj.uorf_start < transcript_obj.uorf_end:
                logging.warning(f"Suspicious negative strand coordinates: uORF_start < uORF_end for {transcript_id}")
            
            # Check for invalid coordinates
            if transcript_obj.uorf_start < 1 or transcript_obj.uorf_end < 1:
                logging.warning(f"Invalid transcript coordinates: less than 1 for {transcript_id}")

    def _reverse_complement(self, sequence: str) -> str:
        """
        Get reverse complement of a DNA sequence.
        
        Args:
            sequence: DNA sequence to complement
            
        Returns:
            Reverse complement of the input sequence
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        """
        Extract transcript ID from BED name field.
        
        Args:
            bed_name: Name field from BED file
            
        Returns:
            Extracted transcript ID
        """
        return bed_name.split('|')[0].split('.')[0]