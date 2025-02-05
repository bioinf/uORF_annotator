import logging
from typing import Optional, Dict
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact


class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile):
        """
        Initialize processor with converter and FASTA reference.
        
        Args:
            converter: CoordinateConverter instance
            fasta: Pysam FASTA file object
        """
        self.converter = converter
        self.fasta = fasta
        self.annotator = None

    def _get_sequence_context(self, chrom: str, start: int, end: int) -> str:
        """
        Get DNA sequence for a genomic region.
        
        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (exclusive)
            
        Returns:
            DNA sequence string
        """
        try:
            return self.fasta.fetch(chrom, start, end).upper()
        except Exception as e:
            logging.error(f"Error fetching sequence for {chrom}:{start}-{end}: {e}")
            return ""

    def process_variant(self, row: pd.Series) -> Optional[Dict]:
        try:
            # Clear previous annotator
            self.annotator = None

            # Extract variant information
            chrom = row['col0']
            vcf_pos = int(row['col1'])
            ref_allele = row['col3']
            alt_allele = row['col4']
            uorf_start = int(row['col9'])
            uorf_end = int(row['col10'])
            strand = row['col13']
            
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
            transcript_pos = self.converter.genome_to_transcript_pos(
                matched_transcript, vcf_pos)

            # Convert 0-based to 1-based position (if position is valid)
            if transcript_pos != "NA":
                transcript_pos = int(transcript_pos) + 1

            # Get sequence and initialize annotator
            sequence = ""
            cds_start = transcript_obj.mainorf_start
            cds_end = transcript_obj.mainorf_end

            # Determine the order of exons based on strand
            exons_to_process = transcript_obj.exons if transcript_obj.strand == '+' else reversed(transcript_obj.exons)
            
            for exon in exons_to_process:
                if exon.genome_end < cds_start or exon.genome_start > cds_end:
                    continue
                
                start = max(exon.genome_start, cds_start)
                end = min(exon.genome_end, cds_end)
                
                exon_seq = self._get_sequence_context(
                    chrom, 
                    start - 1,
                    end
                )
                
                # For negative strand, we need to reverse complement the sequence
                if transcript_obj.strand == '-':
                    exon_seq = self._reverse_complement(exon_seq)
                
                sequence += exon_seq
            self.annotator = VariantAnnotator(sequence, self.fasta)
            
            # Prepare variant data
            variant_data = {
                'position': vcf_pos,
                'ref_allele': ref_allele,
                'alt_allele': alt_allele,
                'uorf_start': uorf_start + 1,
                'uorf_end': uorf_end,
                'strand': strand,
                'maincds_start': transcript_obj.mainorf_start,
                'maincds_end': transcript_obj.mainorf_end,
                'chromosome': chrom
            }

            # Get consequences and impacts
            uorf_consequence = self.annotator.get_consequence(variant_data)
            maincds_impact = None
            if uorf_consequence:
                maincds_impact = self.annotator.predict_impact(variant_data, uorf_consequence)
            
            # Get codon change
            codon_change = self.annotator.get_codon_change(variant_data)
            
            # Check overlap
            overlaps_maincds = self.annotator.does_overlap_maincds(
                uorf_end, transcript_obj.mainorf_start
            )
            
            return {
                'Chromosome': chrom,
                'Original_Genome_Position': vcf_pos,
                'Ref_Allele': ref_allele,
                'Alt_Allele': alt_allele,
                'Strand': strand,
                'Transcript_ID': matched_transcript,
                'Transcript_Position': transcript_pos,
                'uORF_Start': uorf_start,
                'uORF_End': uorf_end,
                'mainCDS_Start': transcript_obj.mainorf_start,
                'mainCDS_End': transcript_obj.mainorf_end,
                'Codon_Change': codon_change,
                'uORF_Consequence': uorf_consequence.value if uorf_consequence else 'None',
                'uORF_mainCDS_Overlap': 'overlapping' if overlaps_maincds else 'non_overlapping',
                'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None'
            }
            
        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}")
            return None

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        """Extract transcript ID from BED name field."""
        return bed_name.split('|')[0].split('.')[0]