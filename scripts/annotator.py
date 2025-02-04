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

    def __init__(self, full_sequence: str):
        self.full_sequence = full_sequence

    def _find_stop_codons(self, sequence: str, frame: int = 0) -> List[int]:
        """Find all stop codons in given frame."""
        stop_positions = []
        for i in range(frame, len(sequence)-2, 3):
            codon = sequence[i:i+3].upper()
            if codon in self.STOP_CODONS:
                stop_positions.append(i)
        return stop_positions

    def _get_sequence_between_features(self, start: int, end: int, strand: str) -> str:
        """Get sequence between two genomic positions."""
        if strand == '+':
            return self.full_sequence[start:end]
        else:
            sequence = self.full_sequence[end:start]
            return self._reverse_complement(sequence)

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))

    def _check_frame_compatibility(self, uorf_end: int, maincds_start: int) -> bool:
        """Check if uORF and mainCDS are in the same frame."""
        return (maincds_start - uorf_end) % 3 == 0

    def get_consequence(self, variant_data: Dict) -> Optional[UORFConsequence]:
        try:
            pos = int(variant_data['position'])
            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']
            uorf_start = int(variant_data['uorf_start'])
            uorf_end = int(variant_data['uorf_end'])
            strand = variant_data['strand']

            if self._affects_start_codon(pos, ref, alt, uorf_start, strand):
                return UORFConsequence.START_LOST

            if self._affects_stop_codon(pos, ref, alt, uorf_end, strand):
                if self._creates_stop_codon(ref, alt):
                    return UORFConsequence.STOP_GAINED
                return UORFConsequence.STOP_LOST

            if self._is_frameshift(ref, alt):
                return UORFConsequence.FRAMESHIFT

            if len(ref) == len(alt) == 1:
                ref_codon = self._get_codon_at_position(pos)
                alt_codon = self._get_altered_codon(ref_codon, pos % 3, ref, alt)
                
                if self._is_synonymous(ref_codon, alt_codon):
                    return UORFConsequence.SYNONYMOUS
                return UORFConsequence.MISSENSE

            if self._affects_splice_region(pos, uorf_start, uorf_end):
                return UORFConsequence.SPLICE_REGION

            return UORFConsequence.NON_CODING

        except Exception as e:
            logging.error(f"Error determining consequence: {str(e)}")
            return None

    def predict_impact(self, variant_data: Dict, uorf_consequence: UORFConsequence) -> MainCDSImpact:
        try:
            uorf_start = int(variant_data['uorf_start'])
            uorf_end = int(variant_data['uorf_end'])
            maincds_start = int(variant_data['maincds_start'])
            maincds_end = int(variant_data['maincds_end'])
            strand = variant_data['strand']
            
            # Simple cases first
            if uorf_consequence in [UORFConsequence.MISSENSE, UORFConsequence.SYNONYMOUS]:
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            overlaps_maincds = self.does_overlap_maincds(uorf_end, maincds_start)

            # Get correct sequence for analysis
            if strand == '+':
                sequence = self._get_sequence_between_features(uorf_end, maincds_end, strand)
                intervening_sequence = self._get_sequence_between_features(uorf_end, maincds_start, strand)
            else:
                sequence = self._get_sequence_between_features(maincds_start, uorf_start, strand)
                intervening_sequence = self._get_sequence_between_features(maincds_end, uorf_end, strand)

            if uorf_consequence == UORFConsequence.STOP_LOST:
                stop_positions = self._find_stop_codons(sequence)
                
                # Check if next stop is mainCDS stop
                if not stop_positions:
                    return MainCDSImpact.N_TERMINAL_EXTENSION
                    
                first_stop = stop_positions[0]
                if overlaps_maincds:
                    if first_stop == (maincds_end - maincds_start):
                        return MainCDSImpact.N_TERMINAL_EXTENSION
                    return MainCDSImpact.OVERLAP_EXTENSION

            elif uorf_consequence == UORFConsequence.FRAMESHIFT:
                frame = (maincds_start - uorf_end) % 3 if strand == '+' else (uorf_start - maincds_end) % 3
                stop_positions = self._find_stop_codons(intervening_sequence, frame)
                
                if not stop_positions:
                    if overlaps_maincds:
                        return MainCDSImpact.OUT_OF_FRAME_OVERLAP
                    return MainCDSImpact.MAIN_CDS_UNAFFECTED
                    
                first_stop = stop_positions[0]
                if overlaps_maincds and maincds_start <= first_stop <= maincds_end:
                    if first_stop != (maincds_end - maincds_start):
                        return MainCDSImpact.UORF_PRODUCT_TRUNCATION

            elif uorf_consequence == UORFConsequence.START_LOST:
                if overlaps_maincds:
                    return MainCDSImpact.OVERLAP_ELIMINATION
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            elif uorf_consequence == UORFConsequence.STOP_GAINED:
                if overlaps_maincds:
                    return MainCDSImpact.OVERLAP_TRUNCATION
                return MainCDSImpact.UORF_PRODUCT_TRUNCATION

            return MainCDSImpact.MAIN_CDS_UNAFFECTED

        except Exception as e:
            logging.error(f"Error predicting mainCDS impact: {str(e)}")
            return None

    def does_overlap_maincds(self, uorf_end: int, maincds_start: int) -> bool:
        """Check if uORF overlaps with mainCDS."""
        return uorf_end >= maincds_start

    def _affects_start_codon(self, pos: int, ref: str, alt: str, uorf_start: int, strand: str) -> bool:
        start_codon_range = range(uorf_start, uorf_start + 3)
        return pos in start_codon_range

    def _affects_stop_codon(self, pos: int, ref: str, alt: str, uorf_end: int, strand: str) -> bool:
        stop_codon_range = range(uorf_end - 2, uorf_end + 1)
        return pos in stop_codon_range

    def _creates_stop_codon(self, ref: str, alt: str) -> bool:
        return alt.upper() in self.STOP_CODONS

    def _is_frameshift(self, ref: str, alt: str) -> bool:
        return abs(len(ref) - len(alt)) % 3 != 0

    def _get_codon_at_position(self, pos: int) -> str:
        codon_start = pos - (pos % 3)
        return self.full_sequence[codon_start:codon_start + 3].upper()

    def _get_altered_codon(self, ref_codon: str, codon_pos: int, ref: str, alt: str) -> str:
        return (ref_codon[:codon_pos] + alt + ref_codon[codon_pos + 1:]).upper()

    def _is_synonymous(self, ref_codon: str, alt_codon: str) -> bool:
        return (self.CODON_TABLE.get(ref_codon) == self.CODON_TABLE.get(alt_codon))

    def _affects_splice_region(self, pos: int, uorf_start: int, uorf_end: int) -> bool:
        splice_ranges = [
            range(uorf_start - 3, uorf_start + 3),
            range(uorf_end - 3, uorf_end + 3)
        ]
        return any(pos in splice_range for splice_range in splice_ranges)