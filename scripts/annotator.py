import logging
from typing import Dict, Optional, Tuple, Callable, List
from enum import Enum

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
        self._initialize_impact_rules()

    def _initialize_impact_rules(self):
        """
        Initialize the rule-based system for impact prediction.
        
        This method sets up the handlers for each consequence type and their order.
        New impact types and rules can be added here in the future.
        """
        # Dictionary mapping consequence types to handler functions
        self.impact_handlers = {
            UORFConsequence.MISSENSE: self._handle_missense_synonymous,
            UORFConsequence.SYNONYMOUS: self._handle_missense_synonymous,
            UORFConsequence.START_LOST: self._handle_start_lost,
            UORFConsequence.STOP_LOST: self._handle_stop_lost,
            UORFConsequence.FRAMESHIFT: self._handle_frameshift,
            UORFConsequence.STOP_GAINED: self._handle_stop_gained,
            UORFConsequence.SPLICE_REGION: self._handle_splice_region,
            UORFConsequence.NON_CODING: self._handle_non_coding
        }

    def does_overlap_maincds(self, uorf_end, maincds_start) -> bool:
        """
        Check if uORF overlaps with mainCDS in transcript coordinates.
        
        Args:
            uorf_end: The end position of the uORF in transcript coordinates
            maincds_start: The start position of the mainCDS in transcript coordinates
            
        Returns:
            True if uORF overlaps with mainCDS, False otherwise
        """
        if uorf_end is None or maincds_start is None:
            return False
            
        # In transcript coordinates, overlap logic is the same regardless of strand
        # uORF overlaps with mainCDS if uORF end >= mainCDS start
        return uorf_end >= maincds_start

    def predict_impact(self, variant_data: Dict, uorf_consequence: UORFConsequence) -> Optional[MainCDSImpact]:
        """
        Predict the impact of a variant on the main CDS.
        
        Args:
            variant_data: Dictionary containing variant information
            uorf_consequence: The consequence type for the variant
            
        Returns:
            The predicted impact on the main CDS, or None if prediction fails
        """
        try:
            # Validate required fields
            required_fields = ['uorf_start', 'uorf_end', 'maincds_start', 'maincds_end', 'position']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    logging.warning(f"Missing required field for impact prediction: {field}")
                    return None

            # Use the appropriate handler for this consequence type
            if uorf_consequence in self.impact_handlers:
                return self.impact_handlers[uorf_consequence](variant_data)
            
            # Default case
            return MainCDSImpact.MAIN_CDS_UNAFFECTED

        except Exception as e:
            logging.error(f"Error predicting mainCDS impact: {str(e)}")
            return None

    def _handle_missense_synonymous(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle missense and synonymous variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_start_lost(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle start-loss variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        if overlaps_maincds:
            return MainCDSImpact.OVERLAP_ELIMINATION
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_stop_lost(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle stop-loss variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        if overlaps_maincds:
            return MainCDSImpact.OVERLAP_EXTENSION
        
        new_stop_pos = self._find_new_stop_codon_position(variant_data, UORFConsequence.STOP_LOST)
        if new_stop_pos is None:
            return MainCDSImpact.MAIN_CDS_UNAFFECTED
            
        # If new stop position is before mainCDS start, it's a truncation
        if new_stop_pos < variant_data['maincds_start']:
            return MainCDSImpact.UORF_PRODUCT_TRUNCATION
        
        # Check if uORF end is adjacent to mainCDS start
        if variant_data['uorf_end'] == variant_data['maincds_start'] - 1:
            return MainCDSImpact.N_TERMINAL_EXTENSION
            
        # Otherwise it's an out-of-frame overlap
        return MainCDSImpact.OUT_OF_FRAME_OVERLAP

    def _handle_frameshift(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle frameshift variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        if overlaps_maincds:
            return MainCDSImpact.OUT_OF_FRAME_OVERLAP
        
        new_stop_pos = self._find_new_stop_codon_position(variant_data, UORFConsequence.FRAMESHIFT)
        if new_stop_pos is None:
            return MainCDSImpact.MAIN_CDS_UNAFFECTED
            
        # If new stop position is before mainCDS start, it's a truncation
        if new_stop_pos < variant_data['maincds_start']:
            return MainCDSImpact.UORF_PRODUCT_TRUNCATION
            
        # If new stop exactly matches mainCDS end, it's an N-terminal extension
        if new_stop_pos == variant_data['maincds_end']:
            return MainCDSImpact.N_TERMINAL_EXTENSION
            
        # Otherwise it's an out-of-frame overlap
        return MainCDSImpact.OUT_OF_FRAME_OVERLAP

    def _handle_stop_gained(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle stop-gain variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        # Currently, stop gain in uORF doesn't affect mainCDS
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_splice_region(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle splice region variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        # Currently, splice region variants in uORF don't affect mainCDS
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_non_coding(self, variant_data: Dict) -> MainCDSImpact:
        """
        Handle non-coding variants.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Main CDS impact classification
        """
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _find_new_stop_codon_position(self, variant_data: Dict, 
                                     uorf_consequence: UORFConsequence) -> Optional[int]:
        """
        Find the position of the new stop codon after a frameshift or stop loss variant.
        
        Args:
            variant_data: Dictionary containing variant information
            uorf_consequence: The consequence type for the variant
            
        Returns:
            Position of the new stop codon in transcript coordinates, or None if no stop codon is found.
        """
        try:
            transcript_pos = variant_data['position']
            uorf_start = variant_data['uorf_start']
            uorf_end = variant_data['uorf_end']
            maincds_end = variant_data.get('maincds_end')
            
            # Determine where to start scanning for a new stop codon
            if uorf_consequence == UORFConsequence.STOP_LOST:
                scan_start_pos = uorf_end  # Start scanning from the end of uORF
            elif uorf_consequence == UORFConsequence.FRAMESHIFT:
                scan_start_pos = transcript_pos  # Start scanning from the variant position
            else:
                return None
            
            # Get sequence from scan start to end of transcript
            sequence = self._get_sequence_from_position(scan_start_pos)
            if not sequence:
                return None
            
            # Calculate the relative position and reading frame
            rel_pos = transcript_pos - uorf_start
            frame = rel_pos % 3
                
            # For frameshift, we need to adjust the frame
            if uorf_consequence == UORFConsequence.FRAMESHIFT:
                # Calculate the frameshift amount (length difference between ref and alt alleles)
                ref_allele = variant_data.get('ref_allele', '')
                alt_allele = variant_data.get('alt_allele', '')
                if not ref_allele or not alt_allele:
                    return None
                    
                # Use actual length difference for frame calculation
                shift_amount = (len(alt_allele) - len(ref_allele)) % 3
                if shift_amount == 0:  # No real frameshift effect
                    return None
                    
                # Adjust the frame based on the shift amount
                frame = (frame + shift_amount) % 3
            
            # Find the next in-frame stop codon
            new_stop_pos = self._find_next_stop_codon(sequence, frame)
            if new_stop_pos is None:
                # If no stop codon found, use mainCDS end as approximation
                if maincds_end is not None:
                    return maincds_end
                return None
                    
            # Convert relative position to transcript position
            return scan_start_pos + new_stop_pos
            
        except Exception as e:
            logging.error(f"Error finding new stop codon: {str(e)}")
            return None
            
    def _get_sequence_from_position(self, start_pos: int) -> str:
        """
        Get the sequence from the given position to the end of the transcript.
        
        Args:
            start_pos: Starting position in transcript coordinates
            
        Returns:
            Sequence string from start_pos to end of transcript
        """
        try:
            # We need to access the full transcript sequence
            full_sequence = self.transcript_seq.sequence
            if not full_sequence:
                return ""
                
            # Calculate the index in the sequence (position is 1-based, index is 0-based)
            start_index = start_pos - 1
            if start_index < 0 or start_index >= len(full_sequence):
                return ""
                
            return full_sequence[start_index:]
            
        except Exception as e:
            logging.error(f"Error getting sequence from position: {str(e)}")
            return ""
            
    def _find_next_stop_codon(self, sequence: str, frame: int) -> Optional[int]:
        """
        Find the position of the next in-frame stop codon in the sequence.
        
        Args:
            sequence: The nucleotide sequence to scan
            frame: The reading frame (0, 1, or 2)
            
        Returns:
            Position of the stop codon relative to the start of the sequence,
            or None if no stop codon is found
        """
        try:
            if not sequence or frame not in (0, 1, 2):
                return None
                
            # Adjust the starting point based on the frame
            start = frame
            
            # Scan the sequence for stop codons in the specified frame
            for i in range(start, len(sequence) - 2, 3):
                codon = sequence[i:i+3].upper()
                if codon in self.STOP_CODONS:
                    return i
                    
            return None
            
        except Exception as e:
            logging.error(f"Error finding next stop codon: {str(e)}")
            return None

    def get_codon_change(self, variant_data: Dict) -> str:
        """
        Get codon change string for the variant.
        
        Args:
            variant_data: Dictionary containing variant information
            
        Returns:
            String representation of the codon change, or "NA" if determination fails
        """
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

            ref_codon = self.transcript_seq.get_codon_at_position(pos)
            if not ref_codon:
                return "NA"

            # Calculate position within codon - use uorf_start as reference
            rel_pos = pos - self.transcript_seq.transcript.uorf_start
            pos_in_codon = rel_pos % 3

            alt_codon = list(ref_codon)
            # Only need to reverse complement for negative strand
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
        """
        Get reverse complement of sequence.
        
        Args:
            sequence: The DNA sequence to complement
            
        Returns:
            The reverse complement of the input sequence
        """
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base.upper(), base) 
                    for base in reversed(sequence))

    def get_consequence(self, variant_data: Dict) -> Optional[UORFConsequence]:
        """
        Determine the consequence of a variant on the uORF.
        
        Args:
            variant_data: Dictionary containing variant information
            
        Returns:
            The consequence type for the variant, or None if determination fails
        """
        try:
            required_fields = ['position', 'ref_allele', 'alt_allele']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    return None

            pos = int(variant_data['position'])
            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']
            
            uorf_start = self.transcript_seq.transcript.uorf_start
            uorf_end = self.transcript_seq.transcript.uorf_end
            
            # Check if position is within the uORF
            if pos < uorf_start or pos > uorf_end:
                return UORFConsequence.NON_CODING

            if self._is_frameshift(ref, alt):
                return UORFConsequence.FRAMESHIFT

            codon_change = self.get_codon_change(variant_data)
            if codon_change == "NA":
                return UORFConsequence.NON_CODING

            ref_codon, alt_codon = codon_change.split('>')

            # Check start codon effect (always at beginning of uORF)
            is_at_start = pos - uorf_start < 3
            if is_at_start and ref_codon in self.START_CODONS:
                if alt_codon not in self.START_CODONS:
                    return UORFConsequence.START_LOST

            # Check stop codon effect (always at end of uORF)
            is_at_end = uorf_end - pos < 3
            if is_at_end:
                if ref_codon in self.STOP_CODONS and alt_codon not in self.STOP_CODONS:
                    return UORFConsequence.STOP_LOST
                if ref_codon not in self.STOP_CODONS and alt_codon in self.STOP_CODONS:
                    return UORFConsequence.STOP_GAINED

            # Check for stop codon anywhere in the uORF
            if alt_codon in self.STOP_CODONS and ref_codon not in self.STOP_CODONS:
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
        """
        Check if variant causes a frameshift.
        
        Args:
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            True if the variant causes a frameshift, False otherwise
        """
        return abs(len(ref) - len(alt)) % 3 != 0

    def _is_synonymous(self, ref_codon: str, alt_codon: str) -> bool:
        """
        Check if amino acid change is synonymous.
        
        Args:
            ref_codon: Reference codon
            alt_codon: Alternate codon
            
        Returns:
            True if the amino acid is unchanged, False otherwise
        """
        return (self.CODON_TABLE.get(ref_codon) == self.CODON_TABLE.get(alt_codon))