from enum import Enum
from typing import Optional, Dict, Tuple, List, Any
import logging

logger = logging.getLogger(__name__)

class uORFConsequence(Enum):
    """Enumeration of possible consequences for upstream Open Reading Frames."""
    
    # Variants with no impact
    NON_CODING = "non_coding_variant"
    INTRONIC = "intronic_variant"

    # Structural uORF impacts
    FRAMESHIFT = "frameshift"
    START_LOST = "start_lost"
    STOP_LOST = "stop_lost"
    STOP_GAINED = "stop_gained"
    
    # Translation impacts
    UORF_TRUNCATION = "uORF_product_truncation"
    UORF_EXTENSION = "uORF_extension"
    
    # Regulatory changes
    REGULATORY_VARIANT = "regulatory_variant"

class VariantAnnotator:
    """
    Annotate variants for their effect on uORFs and main ORFs.
    """
    
    def __init__(self, coordinate_converter):
        """
        Initialize annotator with coordinate converter.
        
        Args:
            coordinate_converter: Instance of CoordinateConverter
        """
        self.converter = coordinate_converter
        
    def _check_frame_impact(self, ref: str, alt: str) -> int:
        """
        Calculate how a variant affects the reading frame.
        
        Args:
            ref: Reference allele
            alt: Alternative allele
            
        Returns:
            Frame shift size (0 if no frame shift)
        """
        return abs(len(alt) - len(ref)) % 3
        
    def _is_stop_codon(self, seq: str) -> bool:
        """
        Check if sequence is a stop codon.
        
        Args:
            seq: Nucleotide sequence to check
        
        Returns:
            True if sequence is a stop codon, False otherwise
        """
        stop_codons = {'TAA', 'TAG', 'TGA'}
        return seq.upper() in stop_codons
        
    def _is_start_codon(self, seq: str) -> bool:
        """
        Check if sequence is a start codon.
        
        Args:
            seq: Nucleotide sequence to check
        
        Returns:
            True if sequence is a start codon, False otherwise
        """
        return seq.upper() == 'ATG'
    
    def _determine_variant_type(self, ref: str, alt: str) -> uORFConsequence:
        """
        Determine variant type based on sequence changes.
        
        Args:
            ref: Reference allele
            alt: Alternative allele
        
        Returns:
            Variant consequence type
        """
        # Length comparison
        if len(ref) == len(alt):
            # Substitution
            if ref != alt:
                return self._classify_substitution(ref, alt)
        elif len(ref) < len(alt):
            # Insertion
            return uORFConsequence.FRAMESHIFT
        elif len(ref) > len(alt):
            # Deletion
            return uORFConsequence.FRAMESHIFT
        
        return uORFConsequence.NON_CODING

    def _classify_substitution(self, ref: str, alt: str) -> uORFConsequence:
        """
        Classify substitution variant.
        
        Args:
            ref: Reference allele
            alt: Alternative allele
        
        Returns:
            Specific variant consequence
        """
        # Check if it creates or removes start/stop codons
        if self._is_start_codon(ref) and not self._is_start_codon(alt):
            return uORFConsequence.START_LOST
        
        if not self._is_stop_codon(ref) and self._is_stop_codon(alt):
            return uORFConsequence.STOP_GAINED
        
        if self._is_stop_codon(ref) and not self._is_stop_codon(alt):
            return uORFConsequence.STOP_LOST
        
        return uORFConsequence.NON_CODING

    def _determine_main_orf_effect(
        self, 
        variant_type: uORFConsequence, 
        ref: str, 
        alt: str, 
        strand: str
    ) -> str:
        """
        Determine the main ORF effect based on variant type and context.
        
        Args:
            variant_type: Variant consequence type
            ref: Reference allele
            alt: Alternative allele
            strand: Transcript strand
        
        Returns:
            Main ORF effect description
        """
        if variant_type == uORFConsequence.FRAMESHIFT:
            # Check start codon context
            if strand == '+':
                return "N-terminal_extension"
            else:
                return "out-of-frame_overlap"
        
        if variant_type == uORFConsequence.STOP_LOST:
            return "overlap_extension"
        
        if variant_type == uORFConsequence.STOP_GAINED:
            return "overlap_truncation"
        
        if variant_type == uORFConsequence.START_LOST:
            return "overlap_elimination"
        
        return "NA"

    def annotate_variant(self, 
                        transcript_id: str,
                        chromosome: str,
                        position: int,
                        ref: str,
                        alt: str,
                        sequence_context: Dict[str, str]) -> Tuple[uORFConsequence, str]:
        """
        Annotate a variant for its effect on uORF.
        
        Args:
            transcript_id: Transcript identifier
            chromosome: Chromosome name
            position: Genomic position
            ref: Reference allele
            alt: Alternative allele
            sequence_context: Dictionary with reference sequences
            
        Returns:
            Tuple of (uORF consequence, additional effect)
        """
        # Check if position is in exon
        if not self.converter.check_position_in_exons(transcript_id, position):
            return uORFConsequence.INTRONIC, "NA"
        
        # Get transcript information
        tx_info = self.converter.get_transcript_info(transcript_id)
        if not tx_info:
            return uORFConsequence.NON_CODING, "NA"
        
        # Determine variant type
        variant_type = self._determine_variant_type(ref, alt)
        
        # Determine ORF effect
        main_orf_effect = self._determine_main_orf_effect(
            variant_type, 
            ref, 
            alt, 
            tx_info['strand']
        )
        
        return variant_type, main_orf_effect
            
    def batch_annotate_variants(self, variants: list) -> list:
        """
        Annotate multiple variants at once.
        
        Args:
            variants: List of variant dictionaries with required fields
            
        Returns:
            List of annotated variants with types and effects
        """
        results = []
        for variant in variants:
            v_type, m_effect = self.annotate_variant(
                variant['transcript_id'],
                variant['chromosome'],
                variant['position'],
                variant['ref'],
                variant['alt'],
                variant['sequence_context']
            )
            
            # Get variant position in transcript
            transcript_details = self.converter.get_variant_transcript_position(
                variant['transcript_id'], 
                variant['position']
            )
            
            results.append({
                **variant,
                'variant_type': v_type.value,
                'main_orf_effect': m_effect,
                'transcript_position': transcript_details['transcript_position'],
                'uorf_coords': transcript_details['uorf_coords'],
                'main_orf_coords': transcript_details['main_orf_coords']
            })
        
        return results