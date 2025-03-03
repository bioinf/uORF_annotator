import logging

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
    
    def _fix_transcript_coordinates(self):
        """
        Fix transcript coordinates to ensure start < end for processing.
        For negative strand, we maintain the biological meaning but ensure
        correct ordering for array indexing.
        """
        if not hasattr(self.transcript, 'uorf_start') or not hasattr(self.transcript, 'uorf_end'):
            return
            
        # Store original coordinates
        self.original_uorf_start = self.transcript.uorf_start
        self.original_uorf_end = self.transcript.uorf_end
        self.original_mainorf_start = self.transcript.mainorf_start
        self.original_mainorf_end = self.transcript.mainorf_end
        
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
        
        Args:
            fasta: FastaFile object for sequence extraction
            
        Returns:
            Complete transcript sequence from 5' to 3'
        """
        transcript_seq = ""
        
        try:
            sorted_exons = self.transcript.exons
            if self.transcript.strand == '-':
                sorted_exons = sorted(self.transcript.exons, 
                                    key=lambda x: x.genome_start, 
                                    reverse=True)
            
            for i, exon in enumerate(sorted_exons):
                try:
                    fetch_start = exon.genome_start - 1
                    fetch_end = exon.genome_end
                    
                    exon_seq = fasta.fetch(self.chromosome, fetch_start, fetch_end).upper()
                    
                    if self.transcript.strand == '-':
                        exon_seq = self._reverse_complement(exon_seq)
                        
                    transcript_seq += exon_seq
                    
                except Exception as e:
                    logging.error(f"Error fetching sequence for exon {i + 1}: {str(e)}")
                    return ""
                    
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
                logging.error(f"Invalid uORF start coordinate: {start}")
                return ""
                
            if end > len(self.sequence):
                logging.error(f"Invalid uORF end coordinate: {end} (transcript length: {len(self.sequence)})")
                return ""
                
            # Validate extraction length
            if end <= start:
                logging.error(f"Invalid uORF length: end ({end}) <= start ({start})")
                return ""
            
            uorf_seq = self.sequence[start:end]
            
            if len(uorf_seq) == 0:
                logging.error("Extracted uORF sequence is empty")
                return ""
                
            return uorf_seq
            
        except Exception as e:
            logging.error(f"Error extracting uORF region: {str(e)}")
            return ""
            
    def get_codon_at_position(self, transcript_pos):
        """
        Get codon at given transcript position.
        
        Args:
            transcript_pos: Position in transcript coordinates
            
        Returns:
            The codon containing the position, or None if retrieval fails
        """
        if not self.uorf_region:
            logging.error("Cannot get codon: uORF region is empty")
            return None
            
        try:
            # Calculate relative position within uORF
            # In transcript coordinates, we always calculate from the start of uORF
            rel_pos = transcript_pos - self.transcript.uorf_start
            
            # Validate relative position
            if rel_pos < 0:
                logging.error(f"Position {transcript_pos} is before uORF start ({self.transcript.uorf_start})")
                return None
                
            codon_start = (rel_pos // 3) * 3
            codon_end = codon_start + 3

            if codon_start < 0:
                logging.error(f"Invalid codon start: {codon_start}")
                return None
                
            if codon_end > len(self.uorf_region):
                logging.error(f"Codon position out of range: {codon_start}")
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