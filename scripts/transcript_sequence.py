import logging

class TranscriptSequence:
    def __init__(self, transcript_obj, fasta, chromosome):
        self.transcript = transcript_obj
        self.chromosome = chromosome

        self.sequence = self._extract_transcript_sequence(fasta)
        if not self.sequence:
            logging.error(f"Failed to extract transcript sequence for {transcript_obj.transcript_id}")
            self.uorf_region = ""
            return

        self.uorf_region = self._extract_uorf_region()
        if not self.uorf_region:
            logging.error(f"Failed to extract uORF region for {transcript_obj.transcript_id}")
        
    def _extract_transcript_sequence(self, fasta):
        """Extract full transcript sequence in 5' to 3' orientation."""
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
        """Extract uORF sequence from transcript sequence."""
        if not self.sequence:
            logging.error("Cannot extract uORF region: transcript sequence is empty")
            return ""
            
        if (self.transcript.uorf_start is None or 
            self.transcript.uorf_end is None):
            logging.error("Cannot extract uORF region: missing coordinates")
            return ""
            
        try:
            if self.transcript.strand == '+':
                start = self.transcript.uorf_start - 1
                end = self.transcript.uorf_end
            else:
                start = min(self.transcript.uorf_start, self.transcript.uorf_end) - 1
                end = max(self.transcript.uorf_start, self.transcript.uorf_end)

            if start < 0:
                logging.error(f"Invalid uORF start coordinate: {start}")
                return ""
                
            if end > len(self.sequence):
                logging.error(f"Invalid uORF end coordinate: {end}")
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
        """Get codon at given transcript position."""
        if not self.uorf_region:
            logging.error("Cannot get codon: uORF region is empty")
            return None
            
        try:
            # Calculate relative position within uORF
            if self.transcript.strand == '+':
                rel_pos = transcript_pos - self.transcript.uorf_start
            else:
                # For negative strand, calculate from the end
                rel_pos = abs(transcript_pos - self.transcript.uorf_end)
            
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
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                     'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base.upper(), base) 
                      for base in reversed(sequence))