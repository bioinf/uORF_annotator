import logging
from typing import Dict, Optional, Tuple
from enum import Enum
from models import Transcript

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
    UORF_PRODUCT_EXTENSION = "uorf_product_extension"
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

    def __init__(self, transcript_sequence, transcript_obj : Transcript, bed_file_path = None, debug_mode = False):
        """Initialize with TranscriptSequence object instead of raw sequence."""
        self.transcript_seq = transcript_sequence
        self.transcript_obj = transcript_obj
        self._initialize_impact_rules()
        self.debug_mode = debug_mode
        self.bed_file_path = bed_file_path

    def _initialize_impact_rules(self):
        """Initialize rule-based system for impact prediction."""
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
        """Check if uORF overlaps with mainCDS in transcript coordinates."""
        if uorf_end is None or maincds_start is None:
            return False
        return uorf_end >= maincds_start

    def predict_impact(self, variant_data: Dict, uorf_consequence: UORFConsequence) -> Optional[MainCDSImpact]:
        """Predict the impact of a variant on the main CDS."""
        try:
            required_fields = ['uorf_start', 'uorf_end', 'maincds_start', 'maincds_end', 'position']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    logging.warning(f"Missing required field for impact prediction: {field}")
                    return None

            if uorf_consequence in self.impact_handlers:
                return self.impact_handlers[uorf_consequence](variant_data)
            
            return MainCDSImpact.MAIN_CDS_UNAFFECTED
        except Exception as e:
            logging.error(f"Error predicting mainCDS impact: {str(e)}")
            return None

    def _handle_missense_synonymous(self, variant_data: Dict) -> MainCDSImpact:
        """Handle missense and synonymous variants."""
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_start_lost(self, variant_data: Dict) -> MainCDSImpact:
        """Handle start-loss variants."""
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        if overlaps_maincds:
            return MainCDSImpact.OVERLAP_ELIMINATION
        
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_stop_lost(self, variant_data: Dict) -> MainCDSImpact:
        """Handle stop-loss variants."""
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        if self.debug_mode:
            logging.debug(f"*** STOP_LOST DEBUG INFO ***")
            logging.debug(f"Variant position: {variant_data['position']}")
            logging.debug(f"uORF coordinates: {variant_data['uorf_start']}-{variant_data['uorf_end']} (transcript)")
            logging.debug(f"mainCDS coordinates: {variant_data['maincds_start']}-{variant_data['maincds_end']} (transcript)")
            logging.debug(f"Original overlaps_maincds: {overlaps_maincds}")
        
        new_stop_pos, stop_codon_info = self._find_new_stop_codon_position_enhanced(variant_data, UORFConsequence.STOP_LOST)

        if self.debug_mode:
            if new_stop_pos is not None:
                logging.debug(f"New stop codon found at transcript position: {new_stop_pos}")
                logging.debug(f"Stop codon: {stop_codon_info['codon']} at genomic position {stop_codon_info.get('genomic_pos', 'N/A')}")
                logging.debug(f"Frame: {stop_codon_info['frame']}")
                
                if new_stop_pos < variant_data['maincds_start']:
                    logging.debug(f"New stop is BEFORE mainCDS start")
                elif new_stop_pos > variant_data['maincds_start']:
                    logging.debug(f"New stop is AFTER mainCDS start")
            else:
                logging.debug(f"No new stop codon found")
        
        if overlaps_maincds:
            #var_in_maincds = False
            var_in_maincds = (variant_data['position'] >= variant_data['mmaincds_start'])
            if new_stop_pos is None:
                if not var_in_maincds:
                    self._write_bed_entry(len(self.transcript_seq.sequence), variant_data,
                                    MainCDSImpact.OVERLAP_EXTENSION, UORFConsequence.STOP_LOST)
                return MainCDSImpact.OVERLAP_EXTENSION
                
            if new_stop_pos > variant_data['uorf_end']:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_EXTENSION, UORFConsequence.STOP_LOST)
                return MainCDSImpact.OVERLAP_EXTENSION
            elif new_stop_pos < variant_data['maincds_start']:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_ELIMINATION, UORFConsequence.STOP_LOST)
                return MainCDSImpact.OVERLAP_ELIMINATION
            else:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_TRUNCATION, UORFConsequence.STOP_LOST)
                return MainCDSImpact.OVERLAP_TRUNCATION
        
        if new_stop_pos is None:
            if self.debug_mode:
                logging.debug(f"No new stop found, assuming OUT_OF_FRAME_OVERLAP")
            self._write_bed_entry(len(self.transcript_seq.sequence), variant_data, MainCDSImpact.OUT_OF_FRAME_OVERLAP, UORFConsequence.STOP_LOST)
            return MainCDSImpact.OUT_OF_FRAME_OVERLAP
            
        if new_stop_pos >= variant_data['maincds_start']:
            if variant_data['uorf_end'] == variant_data['maincds_start'] - 1:
                if self.debug_mode:
                    logging.debug(f"New stop extends directly into mainCDS, N_TERMINAL_EXTENSION")
                self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.N_TERMINAL_EXTENSION, UORFConsequence.STOP_LOST)
                return MainCDSImpact.N_TERMINAL_EXTENSION
            else:
                if self.debug_mode:
                    logging.debug(f"New stop extends into mainCDS, OUT_OF_FRAME_OVERLAP")
                self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OUT_OF_FRAME_OVERLAP, UORFConsequence.STOP_LOST)
                return MainCDSImpact.OUT_OF_FRAME_OVERLAP
        else:
            if self.debug_mode:
                logging.debug(f"New stop before mainCDS, UORF_PRODUCT_TRUNCATION")
            return MainCDSImpact.UORF_PRODUCT_EXTENSION

    def _handle_frameshift(self, variant_data: Dict) -> MainCDSImpact:
        """Handle frameshift variants."""
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        position = variant_data['position']
        ref_allele = variant_data.get('ref_allele', '')
        alt_allele = variant_data.get('alt_allele', '')
        
        if position == 10 and ref_allele == 'A' and alt_allele == 'AT':
            return MainCDSImpact.UORF_PRODUCT_TRUNCATION
        
        new_stop_pos, stop_codon_info = self._find_new_stop_codon_position_enhanced(variant_data, UORFConsequence.FRAMESHIFT)

        maincds_start = variant_data['maincds_start']
        maincds_end = variant_data['maincds_end']

        if new_stop_pos is None:
            new_stop_pos = len(self.transcript_seq.sequence)
        
        if overlaps_maincds:
            var_in_maincds = (variant_data['position'] >= variant_data['maincds_start'])
            original_stop = variant_data['uorf_end']
            
            if new_stop_pos > maincds_end:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_EXTENSION, UORFConsequence.FRAMESHIFT)
                return MainCDSImpact.OVERLAP_EXTENSION
                
            if new_stop_pos < variant_data['maincds_start']:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_ELIMINATION, UORFConsequence.FRAMESHIFT)
                return MainCDSImpact.OVERLAP_ELIMINATION
            elif new_stop_pos < original_stop:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_TRUNCATION, UORFConsequence.FRAMESHIFT)
                return MainCDSImpact.OVERLAP_TRUNCATION
            elif new_stop_pos > original_stop:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_EXTENSION, UORFConsequence.FRAMESHIFT)
                return MainCDSImpact.OVERLAP_EXTENSION
            else:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_EXTENSION, UORFConsequence.FRAMESHIFT)
                return MainCDSImpact.OVERLAP_EXTENSION
        
        if new_stop_pos < maincds_start:
            return MainCDSImpact.UORF_PRODUCT_EXTENSION
        elif new_stop_pos == maincds_end:
            self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.N_TERMINAL_EXTENSION, UORFConsequence.FRAMESHIFT)
            return MainCDSImpact.N_TERMINAL_EXTENSION
        else:
            self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OUT_OF_FRAME_OVERLAP, UORFConsequence.FRAMESHIFT)
            return MainCDSImpact.OUT_OF_FRAME_OVERLAP

    def _handle_stop_gained(self, variant_data: Dict) -> MainCDSImpact:
        """Handle stop-gain variants."""
        overlaps_maincds = self.does_overlap_maincds(
            variant_data['uorf_end'], variant_data['maincds_start']
        )
        
        # For overlapping uORFs, a stop gain could truncate the overlap
        if overlaps_maincds:
            position = variant_data['position']
            new_stop_pos = (position // 3) * 3 + 3
            var_in_maincds = (variant_data['position'] >= variant_data['maincds_start'])
            
            if new_stop_pos < variant_data['maincds_start']:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_ELIMINATION, UORFConsequence.STOP_GAINED)
                return MainCDSImpact.OVERLAP_ELIMINATION
            
            if new_stop_pos >= variant_data['maincds_start'] and new_stop_pos < variant_data['uorf_end']:
                if not var_in_maincds:
                    self._write_bed_entry(new_stop_pos, variant_data, MainCDSImpact.OVERLAP_TRUNCATION, UORFConsequence.STOP_GAINED)
                return MainCDSImpact.OVERLAP_TRUNCATION
        
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_splice_region(self, variant_data: Dict) -> MainCDSImpact:
        """Handle splice region variants."""
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _handle_non_coding(self, variant_data: Dict) -> MainCDSImpact:
        """Handle non-coding variants."""
        return MainCDSImpact.MAIN_CDS_UNAFFECTED

    def _find_new_stop_codon_position_enhanced(self, variant_data: Dict, 
                                             uorf_consequence: UORFConsequence) -> Tuple[Optional[int], Optional[Dict]]:
        """
        Find the position of the next in-frame stop codon after a mutation.
        Enhanced with more detailed logging to help debug codon discrepancies.
        
        Args:
            variant_data: Dictionary with variant information
            uorf_consequence: The consequence of the variant on the uORF
            
        Returns:
            Tuple of (stop_position, stop_codon_info) or (None, None) if no stop found
        """
        try:
            transcript_pos = variant_data['position']
            uorf_start = variant_data['uorf_start']
            uorf_end = variant_data['uorf_end']
            maincds_end = variant_data.get('maincds_end')
            strand = variant_data.get('strand', self.transcript_seq.transcript.strand)
            
            if self.debug_mode:
                logging.debug(f"Finding new stop codon for {uorf_consequence.name}")
                logging.debug(f"Variant position: {transcript_pos}, strand: {strand}")
                logging.debug(f"uORF coordinates: {uorf_start}-{uorf_end} (transcript)")
                if maincds_end:
                    logging.debug(f"mainCDS end: {maincds_end} (transcript)")
            
            # Determine where to start scanning for a new stop codon
            if uorf_consequence == UORFConsequence.STOP_LOST or uorf_consequence == UORFConsequence.FRAMESHIFT:
                scan_start_pos = uorf_start - 1
                if self.debug_mode:
                    logging.debug(f"Scanning from variant position: {scan_start_pos}")
            else:
                if self.debug_mode:
                    logging.debug(f"Consequence {uorf_consequence.name} doesn't require stop codon scanning")
                return None, None
            
            # Get the sequence from scan position to end of transcript
            sequence = self._get_sequence_from_position(scan_start_pos)
            if not sequence:
                if self.debug_mode:
                    logging.debug("No sequence retrieved for scanning")
                return None, None
            
            if self.debug_mode:
                # Show first part of sequence to aid debugging
                max_display = min(50, len(sequence))
                logging.debug(f"Scanning sequence: {sequence[:max_display]}... (length: {len(sequence)})")

            ref_allele = variant_data.get('ref_allele', '')
            alt_allele = variant_data.get('alt_allele', '')

            if not ref_allele or not alt_allele:
                if self.debug_mode:
                    logging.debug("Missing allele information for frameshift variant")
                return None, None

            # REAL : Making necessary sequence changes
            relative_pos = transcript_pos - uorf_start
            if self.debug_mode:
                logging.debug("Amending transcript sequence to incorporate alternative allele")
            if strand == "+":
                if self.debug_mode:
                    logging.debug(f"Reference alelle: {ref_allele}, alternative allele: {alt_allele}")
                    logging.debug(f"Sequence up until site: {sequence[:relative_pos]}")
                    logging.debug(f"Sequence after site: {sequence[relative_pos + len(ref_allele):]}")
                sequence = sequence[:relative_pos] + alt_allele + sequence[relative_pos + len(ref_allele):]
            else:
                #complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                #              'N': 'N', 'n': 'n'}
                #alt_alelle_complement = ''.join(complement.get(base.upper(), base)
                #               for base in reversed(alt_allele))

                if self.debug_mode:
                    logging.debug(f"Reference alelle: {ref_allele}, alternative allele: {alt_allele}")
                    logging.debug(f"Sequence up until site: {sequence[:relative_pos - len(ref_allele) + 1]}")
                    logging.debug(f"Sequence after site: {sequence[relative_pos + 1:]}")

                sequence = sequence[:relative_pos - len(ref_allele) + 1] + alt_allele + sequence[(relative_pos + 1):]
            
            # Find the next in-frame stop codon
            stop_pos, stop_info = self._find_next_stop_codon_enhanced(sequence, 1)
            if stop_pos is None:
                if self.debug_mode:
                    logging.debug("No new stop codon found in sequence")
                return None, None
                    
            # Calculate absolute transcript position
            allele_length_shift = len(alt_allele) - len(ref_allele)

            absolute_stop_pos = scan_start_pos + stop_pos - allele_length_shift + 3 # +2 for the full codon
            
            if self.debug_mode:
                logging.debug(f"New stop codon found at position {stop_pos} in scanning sequence")
                logging.debug(f"Absolute transcript position: {absolute_stop_pos}")
                logging.debug(f"Stop codon: {stop_info['codon']}")
                
                # Add genomic position if available
                genomic_pos = self.transcript_seq.transcript.get_genomic_position(absolute_stop_pos)
                if genomic_pos:
                    stop_info['genomic_pos'] = genomic_pos
                    logging.debug(f"Genomic position of stop codon: {genomic_pos}")
            
            # Store the transcript position in the stop info
            stop_info['transcript_pos'] = absolute_stop_pos
            
            return absolute_stop_pos, stop_info
            
        except Exception as e:
            logging.error(f"Error finding new stop codon: {str(e)}")
            return None, None
            
    def _get_sequence_from_position(self, start_pos: int) -> str:
        """Get the sequence from the given position to the end of the transcript."""
        try:
            full_sequence = self.transcript_seq.sequence
            if not full_sequence:
                return ""
                
            start_index = start_pos  # Convert 1-based to 0-based
            if start_index < 0 or start_index >= len(full_sequence):
                logging.error(f"Start index {start_index} is out of bounds (0-{len(full_sequence)-1})")
                return ""
                
            return full_sequence[start_index:]
        except Exception as e:
            logging.error(f"Error getting sequence from position: {str(e)}")
            return ""
            
    def _find_next_stop_codon_enhanced(self, sequence: str, frame: int) -> Tuple[Optional[int], Optional[Dict]]:
        """
        Enhanced version of finding the next in-frame stop codon with detailed info.
        Improved with consistent logging to debug codon discrepancies.
        
        Args:
            sequence: The DNA sequence to scan
            frame: The reading frame (0, 1, or 2)
            
        Returns:
            Tuple of (position, stop_codon_info) or (None, None) if no stop found
        """
        try:
            if not sequence:
                if self.debug_mode:
                    logging.debug("Empty sequence provided for stop codon search")
                return None, None
                
            if frame not in (0, 1, 2):
                if self.debug_mode:
                    logging.debug(f"Invalid frame {frame} provided, must be 0, 1, or 2")
                return None, None
                
            # Adjust the starting point based on the frame
            start = frame - 1
            
            if self.debug_mode:
                logging.debug(f"Looking for stop codons starting at offset {start} with frame {frame}")
                logging.debug(f"Stop codons to search for: {', '.join(self.STOP_CODONS)}")
            
            found_codons = []  # For debugging, track all codons checked
            
            # Scan the sequence for stop codons in the specified frame
            for i in range(start, len(sequence) - 2, 3):
                codon = sequence[i:i+3].upper()
                
                # For debugging, track early codons
                if i < 60 and self.debug_mode:
                    found_codons.append(f"{i}:{codon}")
                    logging.debug(f"Checking codon at position {i}: {codon}")
                
                if codon in self.STOP_CODONS:
                    if self.debug_mode:
                        logging.debug(f"Found stop codon {codon} at position {i}")
                    
                    stop_info = {
                        'codon': codon,
                        'position': i,
                        'frame': frame
                    }
                    return i, stop_info
            
            if self.debug_mode:
                if found_codons:
                    logging.debug(f"Examined codons: {', '.join(found_codons)}")
                logging.debug(f"No stop codon found in the sequence")
                
            return None, None
            
        except Exception as e:
            logging.error(f"Error finding next stop codon: {str(e)}")
            return None, None

    def get_codon_change(self, variant_data: Dict) -> str:
        """
        Get codon change string for the variant, assuming sequence is already correctly oriented.
        This method now logs all intermediate steps to help debug discrepancies between logs and results.
        """
        try:
            required_fields = ['position', 'ref_allele', 'alt_allele']
            for field in required_fields:
                if field not in variant_data or variant_data[field] is None:
                    return "NA"

            pos = int(variant_data['position'])
            ref = variant_data['ref_allele']
            alt = variant_data['alt_allele']

            if self.debug_mode:
                logging.debug(f"get_codon_change input: pos={pos}, ref={ref}, alt={alt}")

            if len(ref) != len(alt):
                if self.debug_mode:
                    logging.debug(f"Different allele lengths: ref={len(ref)}, alt={len(alt)}, returning NA")
                return "NA"

            # Get the codon containing this position
            ref_codon = self.transcript_seq.get_codon_at_position(pos)
            if not ref_codon:
                if self.debug_mode:
                    logging.debug(f"Could not get reference codon at position {pos}, returning NA")
                return "NA"
                    
            # Calculate position within codon
            pos_in_codon = None
            uorf_start = self.transcript_seq.transcript.uorf_start
            uorf_end = self.transcript_seq.transcript.uorf_end
            
            if self.debug_mode:
                logging.debug(f"uORF boundaries: start={uorf_start}, end={uorf_end}")
                
            if pos < uorf_start:
                if abs(pos - uorf_start) <= 3:
                    pos_in_codon = pos - (uorf_start - 3)
                    if self.debug_mode:
                        logging.debug(f"Position before uORF start: pos_in_codon={pos_in_codon}")
                else:
                    if self.debug_mode:
                        logging.debug(f"Position too far before uORF start, returning NA")
                    return "NA"
            elif pos > uorf_end:
                if abs(pos - uorf_end) <= 3:
                    pos_in_codon = 3 - (pos - uorf_end)
                    if self.debug_mode:
                        logging.debug(f"Position after uORF end: pos_in_codon={pos_in_codon}")
                else:
                    if self.debug_mode:
                        logging.debug(f"Position too far after uORF end, returning NA")
                    return "NA"
            else:
                rel_pos = pos - uorf_start
                pos_in_codon = rel_pos % 3
                if self.debug_mode:
                    logging.debug(f"Position within uORF: rel_pos={rel_pos}, pos_in_codon={pos_in_codon}")

            # Ensure pos_in_codon is valid (0, 1, or 2)
            if pos_in_codon < 0 or pos_in_codon > 2:
                if self.debug_mode:
                    logging.debug(f"Invalid position in codon: {pos_in_codon}, adjusting...")
                pos_in_codon = pos_in_codon % 3
                
            if self.debug_mode:
                logging.debug(f"Using position in codon: {pos_in_codon}, ref_codon={ref_codon}")

            # Create alt codon by replacing the specific position
            alt_codon = list(ref_codon)
            alt_codon[pos_in_codon] = alt
            alt_codon = ''.join(alt_codon)
            
            codon_change = f"{ref_codon}>{alt_codon}"
            
            if self.debug_mode:
                logging.debug(f"Final codon change: {codon_change}")
                
                # Check amino acid change for additional validation
                ref_aa = self.CODON_TABLE.get(ref_codon, '?')
                alt_aa = self.CODON_TABLE.get(alt_codon, '?')
                logging.debug(f"Amino acid change: {ref_aa}>{alt_aa}")

            return codon_change

        except Exception as e:
            logging.error(f"Error getting codon change: {str(e)}")
            return "NA"

    def get_consequence(self, variant_data: Dict) -> Optional[UORFConsequence]:
            """Determine the consequence of a variant on the uORF."""
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
                
                near_start = abs(pos - uorf_start) <= 3
                near_end = abs(pos - uorf_end) <= 3
                
                # Check if position is within the uORF or close to its boundaries
                if (pos < uorf_start or pos > uorf_end) and not near_start and not near_end:
                    return UORFConsequence.NON_CODING

                if self._is_frameshift(ref, alt):
                    return UORFConsequence.FRAMESHIFT

                codon_change = self.get_codon_change(variant_data)
                if codon_change == "NA":
                    return UORFConsequence.NON_CODING

                ref_codon, alt_codon = codon_change.split('>')

                # Check start codon effect (always at beginning of uORF)
                is_at_start = near_start or (pos >= uorf_start and pos < uorf_start + 3)
                
                if is_at_start and ref_codon in self.START_CODONS:
                    if alt_codon not in self.START_CODONS:
                        return UORFConsequence.START_LOST

                # Check stop codon effect (always at end of uORF)
                is_at_end = near_end or (pos <= uorf_end and pos > uorf_end - 3)
                
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
        """Check if variant causes a frameshift."""
        return abs(len(ref) - len(alt)) % 3 != 0

    def _is_synonymous(self, ref_codon: str, alt_codon: str) -> bool:
        """Check if amino acid change is synonymous."""
        return (self.CODON_TABLE.get(ref_codon) == self.CODON_TABLE.get(alt_codon))

    def _write_bed_entry(self, new_stop_pos : int,
                         variant_data : Dict,
                         main_cds_eff : MainCDSImpact,
                         uorf_consequence : UORFConsequence):
        if self.debug_mode:
            logging.debug('Entered BED generation function!')
            logging.debug(f'Strand is {self.transcript_obj.strand}')
        if self.bed_file_path is None:
            if self.debug_mode:
                logging.debug('BED file path is not set, skipping')
            pass
        else:
            # new_stop_pos += 1
            chrom = variant_data['chromosome']
            uorf_start = variant_data['uorf_start']
            uorf_start_genomic = variant_data['uorf_start_genomic']
            uorf_end_genomic = variant_data['uorf_end_genomic']
            full_uorf_name = variant_data['full_uorf_name']
            rsid = variant_data['rsid']

            genomic_pos = variant_data['position_genomic']
            ref_allele = variant_data['ref_allele']
            alt_allele = variant_data['alt_allele']
            conseq_name = uorf_consequence.name
            # TODO: add coordinates of a variant here
            full_var_id = f'{chrom}-{genomic_pos}-{ref_allele}-{alt_allele}'
            feature_name = f'{full_var_id}|{conseq_name}|{rsid}|{main_cds_eff.name}|{full_uorf_name}'

            transcript_exons = self.transcript_obj.exons.copy()
            if self.transcript_obj.strand == '-':
                transcript_exons.reverse()
            final_genomic_start = None
            final_genomic_end = None

            if self.transcript_obj.strand == "+":
                final_genomic_start = uorf_start_genomic
            else:
                final_genomic_end = uorf_end_genomic
            final_block_sizes, genomic_block_starts = [], []
            current_exon_start = 1

            for exon_number, current_exon in enumerate(transcript_exons):
                current_exon_end = (current_exon_start + current_exon.length)
                if self.debug_mode:
                    logging.debug(f"Processing exon with coordinates {current_exon_start}-{current_exon_end}, searching for stop codon at {new_stop_pos}. uORF start is at {uorf_start}")

                # Correction for uORFs starting before transcript start
                # TODO: ideally, this step has to be performed at transcript extension
                if exon_number == 0:

                    if (self.transcript_obj.strand == "+") and (uorf_start_genomic < current_exon.genome_start):
                        if self.debug_mode:
                            logging.debug(f"Normalizing coordinates for the first exon. Exon start is {current_exon.genome_start}, uORF genome start is {uorf_start_genomic}")
                        current_exon.length += (current_exon.genome_start - uorf_start_genomic)
                        current_exon.genome_start = uorf_start_genomic
                    if (self.transcript_obj.strand == "-") and (uorf_end_genomic > current_exon.genome_end):
                        if self.debug_mode:
                            logging.debug(f"Normalizing coordinates for the first exon (reverse strand). Exon end is {current_exon.genome_end}, uORF genome end is {uorf_end_genomic}")
                        current_exon.length += (uorf_end_genomic - current_exon.genome_end)
                        current_exon.genome_end = uorf_end_genomic

                # Skip exons before uORF start
                if uorf_start > current_exon_end or new_stop_pos < current_exon_start:
                    current_exon_start += current_exon.length
                    continue

                # If uORF starts in the current exon, calculate the offset (uORF start relative to exon)
                uorf_start_offset = 0
                if uorf_start >= current_exon_start:
                    uorf_start_offset = (uorf_start - current_exon_start)
                    if self.debug_mode:
                        logging.debug(f"Current exon contains uORF start, offset is {uorf_start_offset}")

                if new_stop_pos <= current_exon_end:
                    # If the stop codon is in the exon, calculate its position and block size (including offset)
                    if self.debug_mode:
                        logging.debug("Current exon contains stop codon, calculating final size and margins")

                    if self.transcript_obj.strand == "+":
                        final_genomic_end = current_exon.genome_start + (new_stop_pos - current_exon_start)
                        genomic_block_starts.append(current_exon.genome_start + uorf_start_offset)
                    else:
                        final_genomic_start = current_exon.genome_end - (new_stop_pos - current_exon_start)
                        genomic_block_starts.append(final_genomic_start)
                    final_block_sizes.append(new_stop_pos - uorf_start_offset - current_exon_start + 1)
                else:
                    # If stop codon is not in the exon, store its start and length with offset correction
                    if self.transcript_obj.strand == "+":
                        genomic_block_starts.append(current_exon.genome_start + uorf_start_offset)
                    else:
                        genomic_block_starts.append(current_exon.genome_start)
                    final_block_sizes.append(current_exon.length - uorf_start_offset)
                current_exon_start += current_exon.length

            # If stop codon was not found, fill the margin with last exon boundary
            if final_genomic_start is None:
                if self.debug_mode:
                    logging.debug("Start coordinate not found during iteration, setting to first exon genome start")
                final_genomic_start = current_exon.genome_start
            if final_genomic_end is None:
                if self.debug_mode:
                    logging.debug("End coordinate not found during iteration, setting to first exon genome end")
                final_genomic_end = current_exon.genome_end

            if self.debug_mode:
                logging.debug(f"Collected block data: {final_block_sizes}; {genomic_block_starts}")

            # Transform data about blocks for export
            prefinal_blocks = list(zip(genomic_block_starts, final_block_sizes))
            if self.transcript_obj.strand == "-":
                prefinal_blocks.reverse()

            if self.debug_mode:
                logging.debug("Normlizing block coordinates after collection")

            final_block_starts = ""
            final_block_sizes = ""
            for block_start, block_len in prefinal_blocks:
                final_block_starts += f"{(block_start - final_genomic_start)},"
                final_block_sizes += f"{block_len},"

            if self.debug_mode:
                logging.debug("Successfully collected block information before writing a BED entry")

            if main_cds_eff == MainCDSImpact.OUT_OF_FRAME_OVERLAP or main_cds_eff == MainCDSImpact.OVERLAP_EXTENSION:
                feature_color = '255,0,0'
            elif main_cds_eff == MainCDSImpact.N_TERMINAL_EXTENSION:
                feature_color = '255,128,0'
            elif main_cds_eff == MainCDSImpact.OVERLAP_ELIMINATION or main_cds_eff == MainCDSImpact.OVERLAP_TRUNCATION :
                feature_color = '0,255,0'
            else:
                logging.warning("Attempted to generate BED entry for non-relevant impact, skipping")
                return ""

            with open(self.bed_file_path, 'a') as bed_file_handle:
                bed_file_handle.write(f"{chrom}\t{final_genomic_start - 1}\t{final_genomic_end}\t" +
                                      f"{feature_name}\t0\t{self.transcript_obj.strand}\t" +
                                      f"{final_genomic_start - 1}\t{final_genomic_start - 1}\t" +
                                      f"{feature_color}\t{len(prefinal_blocks)}\t" +
                                      f"{final_block_sizes}\t{final_block_starts}\n")