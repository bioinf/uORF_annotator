# variant_annotator.py

from enum import Enum
from dataclasses import dataclass
from typing import List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

class VariantType(Enum):
    """Types of variants."""
    SNP = "SNP"
    INSERTION = "INS"
    DELETION = "DEL"
    COMPLEX = "COMPLEX"

class AnnotationType(Enum):
    """Types of variant annotations."""
    STOP_GAINED_CDS_UNAFFECTED = "stop_gained_cds_unaffected"
    FRAMESHIFT_CDS_UNAFFECTED = "frameshift_cds_unaffected"
    SYNONYMOUS_CDS_UNAFFECTED = "synonymous_cds_unaffected"
    MISSENSE_CDS_UNAFFECTED = "missense_cds_unaffected"
    INFRAME_DELETION_CDS_UNAFFECTED = "inframe_deletion_cds_unaffected"
    INFRAME_INSERTION_CDS_UNAFFECTED = "inframe_insertion_cds_unaffected"
    STOP_LOST_CDS_UNAFFECTED = "stop_lost_cds_unaffected"
    STOP_LOST_ORF_OVERLAP = "stop_lost_orf_overlap"
    FRAMESHIFT_ORF_OVERLAP = "frameshift_orf_overlap"
    STOP_LOST_N_TERMINAL_EXTENSION = "stop_lost_n_terminal_extension"
    FRAMESHIFT_N_TERMINAL_EXTENSION = "frameshift_n_terminal_extension"
    NO_STOP_FOUND = "no_stop_found"
    INVALID_REFERENCE = "invalid_reference"
    INTRONIC_VARIANT = "intronic_variant"  # Новый тип для вариантов в интронах

    @classmethod
    def add_annotation_type(cls, name: str, value: str):
        """Add new annotation type dynamically."""
        cls._value2member_map_[value] = cls._member_map_[name] = cls(value)

@dataclass
class ExonBoundary:
    """Representation of exon boundary in transcript coordinates."""
    start: int
    end: int
    examined_bp: int = 0

@dataclass
class Variant:
    """Representation of a genomic variant."""
    ref: str
    alt: str
    pos: int
    strand: str
    
    def get_type(self) -> VariantType:
        """Determine variant type."""
        if len(self.ref) == 1 and len(self.alt) == 1:
            return VariantType.SNP
        elif len(self.ref) < len(self.alt):
            return VariantType.INSERTION
        elif len(self.ref) > len(self.alt):
            return VariantType.DELETION
        else:
            return VariantType.COMPLEX

class VariantAnnotator:
    """Annotates uORF variants based on their effect on protein sequence."""
    
    def __init__(self, transcript_sequence: str, uorf_start: int, uorf_stop: int,
                 cds_start: int, cds_stop: int, exon_boundaries: List[ExonBoundary]):
        self.transcript_sequence = transcript_sequence
        self.uorf_start = uorf_start
        self.uorf_stop = uorf_stop
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.exon_boundaries = exon_boundaries

    def _validate_reference_allele(self, variant: Variant, transcript_pos: int) -> bool:
        """Validate that reference allele matches transcript sequence."""
        try:
            ref_seq = self.transcript_sequence[transcript_pos-1:transcript_pos-1+len(variant.ref)]
            
            # Подробное логирование
            logger.info(f"Reference validation for position {transcript_pos}:")
            logger.info(f"  Full sequence: {self.transcript_sequence}")
            logger.info(f"  Reference allele: {variant.ref}")
            logger.info(f"  Found in sequence: {ref_seq}")
            logger.info(f"  Sequence context: {self.transcript_sequence[max(0, transcript_pos-20):min(len(self.transcript_sequence), transcript_pos+20)]}")
            logger.info(f"  Strand: {variant.strand}")

            return ref_seq.upper() == variant.ref.upper()
        except IndexError:
            logger.warning(f"Position {transcript_pos} is out of bounds for sequence "
                         f"of length {len(self.transcript_sequence)}")
            return False

    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))

    def _find_next_stop_codon(self, sequence: str, start_pos: int) -> Optional[int]:
        """Find next in-frame stop codon in sequence."""
        stop_codons = {'TAG', 'TAA', 'TGA'}
        try:
            logger.debug(f"Searching for stop codon from position {start_pos}")
            logger.debug(f"Sequence from start: {sequence[start_pos:start_pos+30]}")  # Show first 10 codons
            for i in range(start_pos, len(sequence)-2, 3):
                codon = sequence[i:i+3].upper()
                logger.debug(f"Checking codon at position {i}: {codon}")
                if codon in stop_codons:
                    return i
            logger.warning(f"No stop codon found in sequence from position {start_pos}")
            return None
        except IndexError:
            logger.warning(f"Position {start_pos} is out of bounds for sequence "
                         f"of length {len(sequence)}")
            return None

    def _track_examined_bases(self, stop_pos: int):
        """Track number of examined bases in each exon up to stop position."""
        current_pos = 0
        for exon in self.exon_boundaries:
            if current_pos + (exon.end - exon.start) >= stop_pos:
                exon.examined_bp = stop_pos - current_pos
                break
            exon.examined_bp = exon.end - exon.start
            current_pos += exon.end - exon.start

    def _create_modified_sequence(self, variant: Variant, transcript_pos: int) -> str:
        """Create modified transcript sequence with variant."""
        return (self.transcript_sequence[:transcript_pos-1] + 
                variant.alt + 
                self.transcript_sequence[transcript_pos-1+len(variant.ref):])

    def annotate_variant(self, variant: Variant, transcript_pos: int) -> Tuple[AnnotationType, Optional[List[ExonBoundary]]]:
        """Annotate variant based on its effect on uORF and CDS."""
        # Validate reference allele
        if not self._validate_reference_allele(variant, transcript_pos):
            return AnnotationType.INVALID_REFERENCE, None

        # Create modified sequence
        mod_sequence = self._create_modified_sequence(variant, transcript_pos)
        
        # Find new stop codon
        new_stop_pos = self._find_next_stop_codon(mod_sequence, self.uorf_start)
        
        # Детальное логирование для отладки
        logger.info(f"Annotation analysis:")
        logger.info(f"  Original stop pos: {self.uorf_stop}")
        logger.info(f"  New stop pos: {new_stop_pos}")
        logger.info(f"  CDS start: {self.cds_start}")
        logger.info(f"  CDS stop: {self.cds_stop}")
        
        if new_stop_pos is None:
            # Если стоп-кодон не найден, считаем это как потерю стоп-кодона
            # и проверяем, где потенциально может быть новый стоп
            logger.info(f"No stop codon found - checking sequence length: {len(mod_sequence)}")
            if len(mod_sequence) >= self.cds_start:
                return AnnotationType.STOP_LOST_ORF_OVERLAP, self.exon_boundaries
            else:
                return AnnotationType.STOP_LOST_CDS_UNAFFECTED, None
            
        # Track examined bases
        self._track_examined_bases(new_stop_pos)
        
        # Determine variant type
        variant_type = variant.get_type()
        
        # Annotate based on stop codon position
        if new_stop_pos < self.uorf_stop:
            # Case 3
            if variant_type == VariantType.SNP:
                return AnnotationType.STOP_GAINED_CDS_UNAFFECTED, None
            return AnnotationType.FRAMESHIFT_CDS_UNAFFECTED, None
            
        elif new_stop_pos == self.uorf_stop:
            # Case 4
            if variant_type == VariantType.SNP:
                # Compare translated sequences
                orig_prot = self._translate(self.transcript_sequence[self.uorf_start:self.uorf_stop])
                mod_prot = self._translate(mod_sequence[self.uorf_start:self.uorf_stop])
                if orig_prot == mod_prot:
                    return AnnotationType.SYNONYMOUS_CDS_UNAFFECTED, None
                return AnnotationType.MISSENSE_CDS_UNAFFECTED, None
            else:
                if len(variant.alt) < len(variant.ref):
                    return AnnotationType.INFRAME_DELETION_CDS_UNAFFECTED, None
                return AnnotationType.INFRAME_INSERTION_CDS_UNAFFECTED, None
                
        elif new_stop_pos < self.cds_start:
            # Case 5
            if variant_type == VariantType.SNP:
                return AnnotationType.STOP_LOST_CDS_UNAFFECTED, None
            return AnnotationType.FRAMESHIFT_CDS_UNAFFECTED, None
            
        elif new_stop_pos < self.cds_stop or new_stop_pos > self.cds_stop:
            # Case 6
            if variant_type == VariantType.SNP:
                return AnnotationType.STOP_LOST_ORF_OVERLAP, self.exon_boundaries
            return AnnotationType.FRAMESHIFT_ORF_OVERLAP, self.exon_boundaries
            
        elif new_stop_pos == self.cds_stop:
            # Case 7
            if variant_type == VariantType.SNP:
                return AnnotationType.STOP_LOST_N_TERMINAL_EXTENSION, self.exon_boundaries
            return AnnotationType.FRAMESHIFT_N_TERMINAL_EXTENSION, self.exon_boundaries
            
        else:
            raise ValueError("Unexpected stop codon position")

    @staticmethod
    def _translate(sequence: str) -> str:
        """Translate DNA sequence to protein sequence."""
        genetic_code = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
        
        protein = ''
        for i in range(0, len(sequence)-2, 3):
            codon = sequence[i:i+3].upper()
            protein += genetic_code.get(codon, 'X')
        return protein