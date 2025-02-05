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

    def __init__(self, full_sequence: str, fasta_file):
        self.full_sequence = full_sequence
        self.fasta = fasta_file

    def get_codon_change(self, variant_data: Dict) -> str:
        """
        Get codon change string for the variant.
        
        Args:
            variant_data: Dictionary with variant information containing:
                position: genomic position of variant
                ref_allele: reference allele
                alt_allele: alternative allele
                uorf_start: uORF start position
                strand: DNA strand (+ or -)
        
        Returns:
            String in format "XXX>YYY" where XXX is reference codon and YYY is alternative codon,
            or "NA" for indels/errors
        """
        try:
            pos = int(variant_data['position'])
            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']
            strand = variant_data['strand']
            uorf_start = int(variant_data['uorf_start'])
            
            # For indels return NA
            if len(ref) != len(alt):
                return "NA"
                
            # Get position relative to uORF start to determine reading frame
            rel_pos = pos - uorf_start  # Distance from uORF start
            codon_number = rel_pos // 3  # Which codon (0-based)
            pos_in_codon = rel_pos % 3   # Position within codon (0, 1, or 2)
            
            # Calculate codon genomic coordinates
            codon_start = uorf_start + (codon_number * 3)
            
            # Get reference codon sequence
            ref_codon = self.fasta.fetch(
                variant_data['chromosome'],
                codon_start,
                codon_start + 3
            ).upper()
            
            # For reverse strand, reverse complement
            if strand == '-':
                ref_codon = self._reverse_complement(ref_codon)
            
            # Create alternative codon
            alt_codon = list(ref_codon)
            alt_codon[pos_in_codon] = alt
            alt_codon = ''.join(alt_codon)
            
            logging.debug(f"""
                Position: {pos}
                Relative pos: {rel_pos}
                Codon number: {codon_number}
                Position in codon: {pos_in_codon}
                Reference codon: {ref_codon}
                Alternative codon: {alt_codon}
            """)
            
            return f"{ref_codon}>{alt_codon}"
            
        except Exception as e:
            logging.error(f"Error getting codon change: {str(e)}")
            return "NA"

    def _is_synonymous(self, ref_codon: str, alt_codon: str) -> bool:
        """
        Check if codon change is synonymous.
        
        Args:
            ref_codon: Reference codon sequence
            alt_codon: Alternative codon sequence
        
        Returns:
            True if codons code for same amino acid
        """
        return (self.CODON_TABLE.get(ref_codon) == self.CODON_TABLE.get(alt_codon))

    def get_consequence(self, variant_data: Dict) -> Optional[UORFConsequence]:
        """
        Determine the consequence of a variant on the uORF.
        """
        try:
            pos = int(variant_data['position'])
            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']
            uorf_start = int(variant_data['uorf_start'])
            uorf_end = int(variant_data['uorf_end'])
            strand = variant_data['strand']

            # Check frameshift first
            if self._is_frameshift(ref, alt):
                return UORFConsequence.FRAMESHIFT

            # Get codon change
            codon_change = self.get_codon_change(variant_data)
            if codon_change == "NA":
                return UORFConsequence.NON_CODING

            ref_codon, alt_codon = codon_change.split('>')
            logging.debug(f"Analyzing codon change: {ref_codon} -> {alt_codon}")

            # Start codon check (first codon)
            rel_pos = pos - uorf_start
            if rel_pos < 3 and ref_codon in self.START_CODONS:
                if alt_codon not in self.START_CODONS:
                    return UORFConsequence.START_LOST

            # Stop codon check (last codon)
            dist_from_end = uorf_end - pos
            if dist_from_end < 3:
                if ref_codon in self.STOP_CODONS and alt_codon not in self.STOP_CODONS:
                    return UORFConsequence.STOP_LOST
                if ref_codon not in self.STOP_CODONS and alt_codon in self.STOP_CODONS:
                    return UORFConsequence.STOP_GAINED

            # For single nucleotide variants inside uORF
            if len(ref) == len(alt) == 1 and uorf_start <= pos <= uorf_end:
                if self._is_synonymous(ref_codon, alt_codon):
                    return UORFConsequence.SYNONYMOUS
                return UORFConsequence.MISSENSE

            return UORFConsequence.NON_CODING

        except Exception as e:
            logging.error(f"Error determining consequence: {str(e)}")
            return None

    def predict_impact(self, variant_data: Dict, uorf_consequence: UORFConsequence) -> MainCDSImpact:
        """
        Predict the impact of a variant on the main CDS.
        """
        try:
            uorf_start = int(variant_data['uorf_start'])
            uorf_end = int(variant_data['uorf_end'])
            maincds_start = int(variant_data['maincds_start'])
            maincds_end = int(variant_data['maincds_end'])

            # For missense and synonymous variants
            if uorf_consequence in [UORFConsequence.MISSENSE, UORFConsequence.SYNONYMOUS]:
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            overlaps_maincds = self.does_overlap_maincds(uorf_end, maincds_start)

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

    def does_overlap_maincds(self, uorf_end: int, maincds_start: int) -> bool:
        """Check if uORF overlaps with mainCDS."""
        return uorf_end >= maincds_start

    def _is_frameshift(self, ref: str, alt: str) -> bool:
        """Check if variant causes a frameshift."""
        return abs(len(ref) - len(alt)) % 3 != 0

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))