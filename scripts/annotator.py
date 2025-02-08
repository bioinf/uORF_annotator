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
            """
            try:
                # Validate required fields
                required_fields = ['position', 'ref_allele', 'alt_allele', 'uorf_start', 
                                'uorf_end', 'strand', 'chromosome', 'genomic_start', 
                                'genomic_end']
                for field in required_fields:
                    if field not in variant_data or variant_data[field] is None:
                        logging.warning(f"Missing required field: {field}")
                        return "NA"

                pos = int(variant_data['position'])
                ref = variant_data['ref_allele']
                alt = variant_data['alt_allele']
                strand = variant_data['strand']
                uorf_start = int(variant_data['uorf_start'])
                uorf_end = int(variant_data['uorf_end'])
                chromosome = variant_data['chromosome']
                genomic_start = int(variant_data['genomic_start'])
                genomic_end = int(variant_data['genomic_end'])
                
                if len(ref) != len(alt):
                    return "NA"
                    
                if strand == '+':
                    rel_pos = pos - uorf_start
                    codon_number = rel_pos // 3
                    pos_in_codon = rel_pos % 3
                    
                    codon_start = genomic_start + (codon_number * 3)
                    
                    try:
                        ref_codon = self.fasta.fetch(
                            chromosome,
                            codon_start,
                            codon_start + 3
                        ).upper()
                    except Exception as e:
                        logging.error(f"Error fetching reference codon: {str(e)}")
                        return "NA"
                        
                    alt_codon = list(ref_codon)
                    alt_codon[pos_in_codon] = alt
                    alt_codon = ''.join(alt_codon)
                    
                else:  # negative strand
                    rel_pos = pos - uorf_end
                    codon_number = abs(rel_pos) // 3
                    pos_in_codon = abs(rel_pos) % 3
                    
                    codon_end = genomic_end - (codon_number * 3)
                    
                    try:
                        ref_codon = self.fasta.fetch(
                            chromosome,
                            codon_end - 3,
                            codon_end
                        ).upper()
                    except Exception as e:
                        logging.error(f"Error fetching reference codon: {str(e)}")
                        return "NA"

                    # Reverse position in codon for negative strand
                    pos_in_codon = 2 - pos_in_codon
                    
                    # Make change in codon
                    alt_codon = list(ref_codon)
                    alt_codon[pos_in_codon] = alt
                    alt_codon = ''.join(alt_codon)
                    
                    # Reverse both codons
                    ref_codon = self._reverse_complement(ref_codon)
                    alt_codon = self._reverse_complement(alt_codon)

                logging.debug(f"""
                    Strand: {strand}
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
        return (self.CODON_TABLE.get(ref_codon) == self.CODON_TABLE.get(alt_codon))

    def get_consequence(self, variant_data: Dict) -> Optional[UORFConsequence]:
            """
            Determine the consequence of a variant on the uORF.
            """
            try:
                # Validate required fields
                required_fields = ['position', 'ref_allele', 'alt_allele', 'uorf_start', 
                                'uorf_end', 'strand']
                for field in required_fields:
                    if field not in variant_data or variant_data[field] is None:
                        logging.warning(f"Missing required field for consequence analysis: {field}")
                        return None

                pos = int(variant_data['position'])
                ref = variant_data['ref_allele']
                alt = variant_data['alt_allele']
                uorf_start = int(variant_data['uorf_start'])
                uorf_end = int(variant_data['uorf_end'])
                strand = variant_data['strand']

                # Basic position validation
                if pos < min(uorf_start, uorf_end) or pos > max(uorf_start, uorf_end):
                    logging.warning(f"Position {pos} is outside uORF region ({uorf_start}-{uorf_end})")
                    return UORFConsequence.NON_CODING

                # Check frameshift first
                if self._is_frameshift(ref, alt):
                    return UORFConsequence.FRAMESHIFT

                # Get codon change
                codon_change = self.get_codon_change(variant_data)
                if codon_change == "NA":
                    return UORFConsequence.NON_CODING

                ref_codon, alt_codon = codon_change.split('>')
                logging.debug(f"Analyzing codon change: {ref_codon} -> {alt_codon}")

                if strand == '+':
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
                else:
                    # For negative strand
                    # Start codon check (near uorf_end in transcript coordinates)
                    dist_from_start = abs(pos - uorf_end)
                    if dist_from_start < 3 and ref_codon in self.START_CODONS:
                        if alt_codon not in self.START_CODONS:
                            return UORFConsequence.START_LOST

                    # Stop codon check (near uorf_start in transcript coordinates)
                    dist_from_end = abs(pos - uorf_start)
                    if dist_from_end < 3:
                        # Stop codon check
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

    def predict_impact(self, variant_data: Dict, uorf_consequence: UORFConsequence) -> Optional[MainCDSImpact]:
        """Predict the impact of a variant on the main CDS."""
        try:
            # Validate required fields
            required_fields = ['uorf_start', 'uorf_end', 'maincds_start', 'maincds_end', 'strand']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    logging.warning(f"Missing required field for impact prediction: {field}")
                    return None

            uorf_start = int(variant_data['uorf_start'])
            uorf_end = int(variant_data['uorf_end'])
            maincds_start = int(variant_data['maincds_start'])
            maincds_end = int(variant_data['maincds_end'])
            strand = variant_data['strand']

            # For missense and synonymous variants
            if uorf_consequence in [UORFConsequence.MISSENSE, UORFConsequence.SYNONYMOUS]:
                return MainCDSImpact.MAIN_CDS_UNAFFECTED

            overlaps_maincds = self.does_overlap_maincds(uorf_start, uorf_end, maincds_start, maincds_end, strand)

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

    def does_overlap_maincds(self, uorf_start: int, uorf_end: int, maincds_start: int, maincds_end: int, strand: str) -> bool:
        """Check if uORF overlaps with mainCDS in transcript coordinates.
        
        Args:
            uorf_end: End position of uORF in transcript coordinates
            maincds_start: Start position of mainCDS in transcript coordinates
            maincds_end: End position of mainCDS in transcript coordinates
            strand: DNA strand ('+' or '-')
            
        Returns:
            Boolean indicating if overlap exists
        """
        if uorf_end is None or maincds_start is None or maincds_end is None:
            return False
        if strand == '+':
            return uorf_end >= maincds_start
        else:  # negative strand
            return uorf_start >= maincds_end

    def _is_frameshift(self, ref: str, alt: str) -> bool:
        """Check if variant causes a frameshift."""
        return abs(len(ref) - len(alt)) % 3 != 0

    def _reverse_only(self, sequence: str) -> str:
        """Get reverse of sequence without complementing."""
        return sequence[::-1]

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))