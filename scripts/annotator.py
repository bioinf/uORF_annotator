from enum import Enum
from typing import Optional, Dict, List, Tuple
import logging


class UORFConsequence(Enum):
    START_LOST = "uorf_start_lost"
    STOP_LOST = "uorf_stop_lost"
    STOP_GAINED = "uorf_stop_gained"
    FRAMESHIFT = "uorf_frameshift"
    MISSENSE = "uorf_missense"
    SYNONYMOUS = "uorf_synonymous"
    SPLICE_REGION = "uorf_splice_region"
    NON_CODING = "non_coding"


class MainCDSImpact(Enum):
    N_TERMINAL_EXTENSION = "n_terminal_extension"
    OUT_OF_FRAME_OVERLAP = "out_of_frame_overlap"
    UORF_PRODUCT_TRUNCATION = "uorf_product_truncation"
    STOP_GAINED = "stop_gained"
    OVERLAP_EXTENSION = "overlap_extension"
    OVERLAP_TRUNCATION = "overlap_truncation"
    OVERLAP_ELIMINATION = "overlap_elimination"
    MAIN_CDS_UNAFFECTED = "main_cds_unaffected"


class VariantAnnotator:
    CODON_TABLE = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    
    START_CODONS = {'ATG'}
    STOP_CODONS = {'TAA', 'TAG', 'TGA'}

    def __init__(self, transcript_sequence):
        """
        Initialize with TranscriptSequence object instead of raw sequence and fasta.
        
        Args:
            transcript_sequence (TranscriptSequence): Object containing transcript sequence data
        """
        self.transcript_seq = transcript_sequence

    def does_overlap_maincds(self, transcript_strand: str, uorf_start_genomic: int, uorf_end_genomic: int, 
                        maincds_start_genomic: int, maincds_end_genomic: int) -> bool:
        """
        Check if uORF overlaps with mainCDS in genomic coordinates.
        For + strand, checks if uORF end is after mainCDS start.
        For - strand, checks if uORF start is before mainCDS end.
        """
        if uorf_start_genomic is None or uorf_end_genomic is None or maincds_start_genomic is None or maincds_end_genomic is None:
            return False
            
        if transcript_strand == '+':
            return uorf_end_genomic >= maincds_start_genomic
        else:  # strand == '-'
            return uorf_start_genomic <= maincds_end_genomic

    def predict_impact(self, variant_data: Dict, uorf_consequence: UORFConsequence) -> Optional[MainCDSImpact]:
        """Predict the impact of a variant on the main CDS."""
        try:
            # For missense and synonymous variants
            if uorf_consequence in [UORFConsequence.MISSENSE, UORFConsequence.SYNONYMOUS]:
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            required_fields = ['overlaps_maincds']
            for field in required_fields:
                if field not in variant_data:
                    logging.warning(f"Missing required field for impact prediction: {field}")
                    return None

            overlaps_maincds = variant_data['overlaps_maincds']

            # For start loss
            if uorf_consequence == UORFConsequence.START_LOST:
                if overlaps_maincds:
                    return MainCDSImpact.OVERLAP_ELIMINATION
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            # For stop loss
            if uorf_consequence == UORFConsequence.STOP_LOST:
                if overlaps_maincds:
                    return MainCDSImpact.OVERLAP_EXTENSION
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            # For frameshift
            if uorf_consequence == UORFConsequence.FRAMESHIFT:
                if overlaps_maincds:
                    return MainCDSImpact.OUT_OF_FRAME_OVERLAP
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            return MainCDSImpact.MAIN_CDS_UNAFFECTED

        except Exception as e:
            logging.error(f"Error predicting mainCDS impact: {str(e)}")
            return None

    def get_codon_change(self, variant_data: Dict) -> str:
        """Get codon change string for the variant."""
        try:
            required_fields = ['position', 'ref_allele', 'alt_allele']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    return "NA"

            pos = int(variant_data['position'])
            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']

            if len(ref) != len(alt):
                return "NA"

            # Get reference codon using TranscriptSequence logic
            ref_codon = self.transcript_seq.get_codon_at_position(pos)
            if not ref_codon:
                return "NA"
                
            # Calculate relative position using same logic as in TranscriptSequence
            if self.transcript_seq.transcript.strand == '+':
                rel_pos = pos - self.transcript_seq.transcript.uorf_start
            else:
                rel_pos = abs(pos - self.transcript_seq.transcript.uorf_end)
                
            # Calculate position within codon
            pos_in_codon = rel_pos % 3

            # Create alternate codon
            alt_codon = list(ref_codon)
            if self.transcript_seq.transcript.strand == '-':
                alt = self._reverse_complement(alt)
            alt_codon[pos_in_codon] = alt
            alt_codon = ''.join(alt_codon)

            return f"{ref_codon}>{alt_codon}"

        except Exception as e:
            logging.error(f"Error getting codon change: {str(e)}")
            return "NA"


    @staticmethod
    def _reverse_complement(sequence):
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base.upper(), base) 
                    for base in reversed(sequence))

    def get_consequence(self, variant_data: Dict) -> Optional[UORFConsequence]:
        """Determine the consequence of a variant on the uORF."""
        try:
            required_fields = ['position', 'ref_allele', 'alt_allele', 'codon_change']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    return None

            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']

            if self._is_frameshift(ref, alt):
                return UORFConsequence.FRAMESHIFT

            # Use pre-calculated codon change
            codon_change = variant_data['codon_change']
            if codon_change == "NA":
                return UORFConsequence.NON_CODING

            ref_codon, alt_codon = codon_change.split('>')

            # Check start codon effect
            if ref_codon in self.START_CODONS:
                if alt_codon not in self.START_CODONS:
                    return UORFConsequence.START_LOST

            # Check stop codon effect
            if ref_codon in self.STOP_CODONS and alt_codon not in self.STOP_CODONS:
                return UORFConsequence.STOP_LOST
            if ref_codon not in self.STOP_CODONS and alt_codon in self.STOP_CODONS:
                return UORFConsequence.STOP_GAINED

            # For single nucleotide variants inside uORF
            if len(ref) == len(alt) == 1:
                if self._is_synonymous(ref_codon, alt_codon):
                    return UORFConsequence.SYNONYMOUS
                return UORFConsequence.MISSENSE

            return UORFConsequence.NON_CODING

        except Exception as e:
            logging.error(f"Error determining consequence: {str(e)}")
            return None

    def _is_frameshift(self, ref: str, alt: str) -> bool:
        """Check if variant causes a frameshift."""
        return abs(len(ref) - len(alt)) % 3 != 0

    def _is_synonymous(self, ref_codon: str, alt_codon: str) -> bool:
        """Check if amino acid change is synonymous."""
        return (self.CODON_TABLE.get(ref_codon) == self.CODON_TABLE.get(alt_codon))
