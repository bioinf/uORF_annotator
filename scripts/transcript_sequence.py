"""
Module for handling transcript sequences and uORF extraction.
"""
import logging
from typing import Optional

class TranscriptSequence:
    def __init__(self, transcript_obj, fasta, chromosome):
        """
        Initialize with transcript object, fasta reference, and chromosome.
        
        Args:
            transcript_obj: Transcript object containing coordinate information
            fasta: FastaFile object containing reference sequence
            chromosome: Chromosome name/identifier
        """
        self.transcript = transcript_obj
        self.chromosome = chromosome

        # Store original coordinates before any modifications
        self._store_original_coordinates()

        # Fix coordinate ordering for negative strand transcripts before extraction
        self._fix_transcript_coordinates()

        self.sequence = self._extract_transcript_sequence(fasta)
        if not self.sequence:
            logging.error(f"Failed to extract transcript sequence for {transcript_obj.transcript_id}")
            self.uorf_region = ""
            return

        self.uorf_region = self._extract_uorf_region()
        if not self.uorf_region:
            logging.error(f"Failed to extract uORF region for {transcript_obj.transcript_id}")
    
    def _store_original_coordinates(self):
        """
        Store the original biological coordinates before any manipulation.
        This ensures we can always refer back to the true biological coordinates.
        """
        if not hasattr(self.transcript, 'uorf_start') or not hasattr(self.transcript, 'uorf_end'):
            return
            
        # Store original coordinates
        self.original_uorf_start = self.transcript.uorf_start
        self.original_uorf_end = self.transcript.uorf_end
        self.original_mainorf_start = self.transcript.mainorf_start
        self.original_mainorf_end = self.transcript.mainorf_end
        
        # Store original strand information
        self.original_strand = self.transcript.strand
    
    def _fix_transcript_coordinates(self):
        """
        Fix transcript coordinates to ensure start < end for processing.
        For negative strand, we maintain the biological meaning but ensure
        correct ordering for array indexing.
        """
        if not hasattr(self.transcript, 'uorf_start') or not hasattr(self.transcript, 'uorf_end'):
            return
        
        # For negative strand, ensure start < end for processing
        if self.transcript.strand == '-':
            # Fix uORF coordinates if needed
            if (self.transcript.uorf_start is not None and 
                self.transcript.uorf_end is not None and 
                self.transcript.uorf_start > self.transcript.uorf_end):
                self.transcript.uorf_start, self.transcript.uorf_end = (
                    self.transcript.uorf_end, self.transcript.uorf_start
                )
                logging.info(f"Fixed uORF coordinates for negative strand transcript {self.transcript.transcript_id}")
            
            # Fix mainORF coordinates if needed
            if (self.transcript.mainorf_start is not None and 
                self.transcript.mainorf_end is not None and 
                self.transcript.mainorf_start > self.transcript.mainorf_end):
                self.transcript.mainorf_start, self.transcript.mainorf_end = (
                    self.transcript.mainorf_end, self.transcript.mainorf_start
                )
                logging.info(f"Fixed mainORF coordinates for negative strand transcript {self.transcript.transcript_id}")
        
        # Ensure coordinates are valid (1-based)
        if self.transcript.uorf_start is not None and self.transcript.uorf_start < 1:
            logging.warning(f"Invalid uORF start position {self.transcript.uorf_start} for {self.transcript.transcript_id}, setting to 1")
            self.transcript.uorf_start = 1
        
        if self.transcript.mainorf_start is not None and self.transcript.mainorf_start < 1:
            logging.warning(f"Invalid mainORF start position {self.transcript.mainorf_start} for {self.transcript.transcript_id}, setting to 1")
            self.transcript.mainorf_start = 1

    def _extract_transcript_sequence(self, fasta):
        """
        Extract full transcript sequence in 5' to 3' orientation.
        Also handle cases where transcript has been extended to include uORF.
        
        Args:
            fasta: FastaFile object for sequence extraction
            
        Returns:
            Complete transcript sequence from 5' to 3'
        """
        transcript_seq = ""
        
        try:
            was_extended = getattr(self.transcript, 'was_extended', False)
            
            # Use extended boundaries for extended transcripts
            if was_extended:
                logging.debug(f"Using extended boundaries for transcript {self.transcript.transcript_id}")
                
                # For extended transcripts, use the entire range between min/max genomic positions
                min_genomic_pos = min(self.transcript.genome_to_transcript.keys())
                max_genomic_pos = max(self.transcript.genome_to_transcript.keys())
                
                # Extract the entire region for extended transcript
                fetch_start = min_genomic_pos - 1  # Convert to 0-based
                fetch_end = max_genomic_pos
                
                try:
                    # Get sequence in genomic coordinates
                    full_seq = fasta.fetch(self.chromosome, fetch_start, fetch_end).upper()
                    
                    # For negative strand, do reverse complement
                    if self.transcript.strand == '-':
                        logging.debug(f"Doing reverse complement for negative strand extended transcript "
                                    f"{self.transcript.transcript_id} (length: {len(full_seq)})")
                        full_seq = self._reverse_complement(full_seq)
                        
                    return full_seq
                except Exception as e:
                    logging.error(f"Error fetching sequence: {str(e)}")
                    return ""
            
            # Regular processing for non-extended transcripts
            sorted_exons = self.transcript.exons
            
            # For negative strand, exons need to be processed in reverse order
            if self.transcript.strand == '-':
                sorted_exons = sorted(self.transcript.exons, 
                                    key=lambda x: x.genome_start, 
                                    reverse=True)
                logging.debug(f"Sorted {len(sorted_exons)} exons for negative strand transcript "
                            f"{self.transcript.transcript_id}")
            
            for i, exon in enumerate(sorted_exons):
                try:
                    fetch_start = exon.genome_start - 1
                    fetch_end = exon.genome_end
                    
                    exon_seq = fasta.fetch(self.chromosome, fetch_start, fetch_end).upper()
                    
                    if self.transcript.strand == '-':
                        exon_seq = self._reverse_complement(exon_seq)
                        logging.debug(f"Reverse complemented exon {i} for negative strand "
                                    f"(length: {len(exon_seq)})")
                        
                    transcript_seq += exon_seq
                    
                except Exception as e:
                    logging.error(f"Error fetching sequence for exon: {str(e)}")
                    return ""
                    
            # Log sequence details
            logging.debug(f"Extracted transcript sequence for {self.transcript.transcript_id}, "
                        f"strand: {self.transcript.strand}, length: {len(transcript_seq)}")
                
            return transcript_seq
            
        except Exception as e:
            logging.error(f"Error in transcript sequence extraction: {str(e)}")
            return ""

    def _extract_uorf_region(self):
        """
        Extract uORF sequence from transcript sequence.
        
        Returns:
            Sequence of the uORF region
        """
        if not self.sequence:
            logging.error("Cannot extract uORF region: transcript sequence is empty")
            return ""
            
        if (self.transcript.uorf_start is None or 
            self.transcript.uorf_end is None):
            logging.error("Cannot extract uORF region: missing coordinates")
            return ""
            
        try:
            # Convert from 1-based to 0-based coordinates for sequence indexing
            start = self.transcript.uorf_start - 1
            end = self.transcript.uorf_end
            
            # Validate coordinates
            if start < 0:
                logging.warning(f"Invalid uORF start coordinate: {start+1} (adjusted to 1)")
                start = 0
                
            if end > len(self.sequence):
                logging.warning(f"Invalid uORF end coordinate: {end} (adjusted to transcript length: {len(self.sequence)})")
                end = len(self.sequence)
                
            # Check if start is still less than end after adjustments
            if end <= start:
                logging.error(f"Invalid uORF length after adjustment: end ({end}) <= start ({start})")
                return ""
            
            uorf_seq = self.sequence[start:end]
            
            # Log details about uORF extraction
            logging.debug(f"Extracted uORF region for {self.transcript.transcript_id}, "
                        f"strand: {self.transcript.strand}, coordinates: {start+1}-{end}, "
                        f"length: {len(uorf_seq)}")
            
            if len(uorf_seq) == 0:
                logging.error("Extracted uORF sequence is empty")
                return ""
                
            # Check if uORF sequence has proper start/stop codons for debugging
            if len(uorf_seq) >= 3:
                start_codon = uorf_seq[:3]
                stop_codon = uorf_seq[-3:] if len(uorf_seq) >= 6 else None
                logging.debug(f"uORF start codon: {start_codon}, stop codon: {stop_codon}")
                
                # Note: Since we've already done reverse complement for negative strand,
                # the start codon should be ATG in the extracted sequence for both strands
                if start_codon != 'ATG':
                    logging.warning(f"Unusual start codon for uORF: {start_codon}, expected 'ATG'")
                    
            return uorf_seq
            
        except Exception as e:
            logging.error(f"Error extracting uORF region: {str(e)}")
            return ""

    def get_codon_at_position(self, transcript_pos):
        """
        Get codon at given transcript position, with correct handling for negative strand and extended transcripts.
        
        Args:
            transcript_pos: Position in transcript coordinates
            
        Returns:
            The codon containing the position, or None if retrieval fails
        """
        if not self.uorf_region:
            logging.error("Cannot get codon: uORF region is empty")
            return None
            
        try:
            was_extended = getattr(self.transcript, 'was_extended', False)
            
            logging.debug(f"Getting codon at position {transcript_pos} for transcript {self.transcript.transcript_id}, "
                        f"extended: {was_extended}")
            
            # Check if position is outside uORF boundaries
            near_start = abs(transcript_pos - self.transcript.uorf_start) <= 3
            near_end = abs(transcript_pos - self.transcript.uorf_end) <= 3
            
            if transcript_pos < self.transcript.uorf_start and near_start:
                # Position is before uORF start but within 3 bp
                codon = self.uorf_region[:3]
                logging.debug(f"Position is near uORF start, using first codon: {codon}")
                return codon
                
            if transcript_pos > self.transcript.uorf_end and near_end:
                # Position is after uORF end but within 3 bp
                codon = self.uorf_region[-3:]
                logging.debug(f"Position is near uORF end, using last codon: {codon}")
                return codon
                
            # For position within uORF, calculate codon based on reading frame
            rel_pos = transcript_pos - self.transcript.uorf_start
            if rel_pos < 0 or rel_pos >= len(self.uorf_region):
                logging.error(f"Position {transcript_pos} is outside uORF range ({self.transcript.uorf_start}-{self.transcript.uorf_end})")
                return None
                
            # Calculate the frame for this position (0, 1, or 2)
            frame = rel_pos % 3
            
            # Calculate the codon start by finding the start of the codon containing this position
            codon_start = rel_pos - frame
            codon_end = codon_start + 3

            # Check if the calculated position is valid
            if codon_start < 0:
                logging.error(f"Invalid codon start: position too close to uORF start")
                return None
                
            if codon_end > len(self.uorf_region):
                logging.error(f"Codon position out of range: codon extends beyond uORF region")
                return None
            
            codon = self.uorf_region[codon_start:codon_end]
            return codon
                
        except Exception as e:
            logging.error(f"Error getting codon: {str(e)}")
            return None

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