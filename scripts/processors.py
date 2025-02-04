# processors.py

import logging
from typing import Optional, Dict
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact


class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta_file: pysam.FastaFile):
        """
        Initialize processor with converter and FASTA reference.
        """
        self.converter = converter
        self.fasta = fasta_file
        self.annotator = None

    def _get_full_sequence_context(self, chrom: str, variant_data: Dict) -> str:
        """Get full sequence context for variant analysis."""
        uorf_start = int(variant_data['uorf_start'])
        uorf_end = int(variant_data['uorf_end'])
        maincds_start = int(variant_data['maincds_start'])
        maincds_end = int(variant_data['maincds_end'])
        strand = variant_data['strand']
        
        if strand == '+':
            seq_start = min(uorf_start, maincds_start)
            seq_end = max(uorf_end, maincds_end)
        else:
            seq_start = min(uorf_start, maincds_start)
            seq_end = max(uorf_end, maincds_end)
                
        return self.fasta.fetch(chrom, seq_start - 3, seq_end + 3)

    def process_variant(self, row: pd.Series) -> Optional[Dict]:
        try:
            # Extract basic variant information
            chrom = row['col0']
            vcf_pos = int(row['col1'])
            ref_allele = row['col3']
            alt_allele = row['col4']
            
            # Get transcript information
            transcript_id = self._extract_transcript_id(row['col11'])
            matching_transcripts = [
                tid for tid in self.converter.transcripts.keys() 
                if tid.startswith(transcript_id)
            ]

            if not matching_transcripts:
                logging.warning(f"Transcript {transcript_id} not found")

                return None

            matched_transcript = matching_transcripts[0]
            transcript_obj = self.converter.transcripts[matched_transcript]
            transcript_pos = self.converter.genome_to_transcript_pos(
                matched_transcript, vcf_pos)

            # Prepare variant data
            variant_data = {
                'position': vcf_pos,
                'ref_allele': ref_allele,
                'alt_allele': alt_allele,
                'uorf_start': row['col9'],
                'uorf_end': row['col10'],
                'strand': row['col13'],
                'maincds_start': transcript_obj.mainorf_start,
                'maincds_end': transcript_obj.mainorf_end
            }

            # Get sequence context and initialize annotator
            full_sequence = self._get_full_sequence_context(chrom, variant_data)
            if not self.annotator:
                self.annotator = VariantAnnotator(full_sequence)

            # Get variant consequences and impacts
            uorf_consequence = self.annotator.get_consequence(variant_data)
            maincds_impact = None
            if uorf_consequence:
                maincds_impact = self.annotator.predict_impact(variant_data, uorf_consequence)

            # Check for overlap with mainCDS
            overlaps_maincds = self.annotator.does_overlap_maincds(
                int(row['col10']), transcript_obj.mainorf_start
            )

            # Return complete variant information
            return {
                'Chromosome': chrom,
                'Original_Genome_Position': vcf_pos,
                'Ref_Allele': ref_allele,
                'Alt_Allele': alt_allele,
                'Strand': row['col13'],
                'Transcript_ID': matched_transcript,
                'Transcript_Position': transcript_pos,
                'uORF_Start': row['col9'],
                'uORF_End': row['col10'],
                'mainCDS_Start': transcript_obj.mainorf_start,
                'mainCDS_End': transcript_obj.mainorf_end,
                'Overlaps_mainCDS': overlaps_maincds,
                'uORF_Consequence': uorf_consequence.value if uorf_consequence else 'None',
                'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None'
            }

        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}")
            return None

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        """Extract transcript ID from BED name field."""
        return bed_name.split('|')[0].split('.')[0]