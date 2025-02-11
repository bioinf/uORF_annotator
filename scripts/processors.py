import logging
from typing import Optional, Dict
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact
from transcript_sequence import TranscriptSequence

class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile):
        self.converter = converter
        self.fasta = fasta
        self.annotator = None

    def process_variant(self, row: pd.Series) -> Optional[Dict]:
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
                    
            # Initialize annotator with only TranscriptSequence
            self.annotator = VariantAnnotator(transcript_seq)

            # Check overlap using genomic coordinates
            overlaps_maincds = (
                self.annotator.does_overlap_maincds(
                    transcript_obj.strand,
                    transcript_obj.uorf_start_genomic,
                    transcript_obj.uorf_end_genomic,
                    transcript_obj.mainorf_start_genomic,
                    transcript_obj.mainorf_end_genomic
                )
                if (transcript_obj.uorf_start_genomic is not None and
                    transcript_obj.uorf_end_genomic is not None and
                    transcript_obj.mainorf_start_genomic is not None and
                    transcript_obj.mainorf_end_genomic is not None)
                else False
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
                'uORF_Start_Transcript': transcript_obj.uorf_start,
                'uORF_End_Transcript': transcript_obj.uorf_end,
                'mainCDS_Start_Genomic': transcript_obj.mainorf_start_genomic,
                'mainCDS_End_Genomic': transcript_obj.mainorf_end_genomic,
                'mainCDS_Start_Transcript': transcript_obj.mainorf_start,
                'mainCDS_End_Transcript': transcript_obj.mainorf_end,
                'Codon_Change': codon_change,
                'uORF_Consequence': uorf_consequence.value if uorf_consequence else 'None',
                'uORF_mainCDS_Overlap': 'overlapping' if overlaps_maincds else 'non_overlapping',
                'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None'
            }
                    
        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}")
            return None

    def _reverse_complement(self, sequence: str) -> str:
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        return bed_name.split('|')[0].split('.')[0]