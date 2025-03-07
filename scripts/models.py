"""
Models for representing genomic and transcript data structures.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, NamedTuple
import logging


class Coordinates(NamedTuple):
    """Store both genomic and transcript coordinates."""
    genomic: int
    transcript: int


@dataclass
class Exon:
    start: int
    length: int
    genome_start: int
    genome_end: int


class Transcript:
    def __init__(self, 
                 transcript_id: str,
                 chromosome: str,
                 strand: str,
                 exons: List[Exon],
                 mainorf_start: Optional[int] = None, 
                 mainorf_end: Optional[int] = None,
                 uorf_start: Optional[int] = None,
                 uorf_end: Optional[int] = None):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = sorted(exons, key=lambda x: x.genome_start)
        
        # Store genomic coordinates
        self.mainorf_start_genomic = mainorf_start
        self.mainorf_end_genomic = mainorf_end
        self.uorf_start_genomic = uorf_start
        self.uorf_end_genomic = uorf_end
        
        # Initialize coordinate maps
        self.genome_to_transcript = {}
        self.transcript_to_genome = {}
        
        # Flag to indicate if the transcript was extended to include uORF
        self.was_extended = False
        
        # Build initial coordinate maps from exons
        self._build_coordinate_maps()
        
        # First check if uORF is outside of transcript boundaries (in genomic coordinates)
        # and extend the transcript if necessary
        self._check_and_extend_for_uorf()
        
        # AFTER extension, convert genomic coordinates to transcript coordinates
        self.mainorf_start = self._convert_to_transcript_coords(mainorf_start)
        self.mainorf_end = self._convert_to_transcript_coords(mainorf_end)
        self.uorf_start = self._convert_to_transcript_coords(uorf_start)
        self.uorf_end = self._convert_to_transcript_coords(uorf_end)

    def _convert_to_transcript_coords(self, genomic_pos: int) -> Optional[int]:
        """Convert genomic position to transcript coordinates."""
        if genomic_pos is None:
            return None
        return self.genome_to_transcript.get(genomic_pos)

    def _build_coordinate_maps(self):
        """Build mappings between genomic and transcript coordinates."""
        if self.strand == '+':
            transcript_pos = 1  # Start counting from 1
            
            for exon in self.exons:
                offset = 0  # Keep track of position within current exon
                for genome_pos in range(exon.genome_start, exon.genome_end + 1):
                    # For first position in exon, use transcript_pos directly
                    current_pos = transcript_pos + offset
                    self.genome_to_transcript[genome_pos] = current_pos
                    self.transcript_to_genome[current_pos] = genome_pos
                    offset += 1
                transcript_pos += offset
        else:
            # For negative strand, start from highest genomic coordinate
            transcript_pos = 1
            sorted_exons = sorted(self.exons, key=lambda x: x.genome_end, reverse=True)
            
            for exon in sorted_exons:
                offset = 0
                # Count positions from end to start for negative strand
                for genome_pos in range(exon.genome_end, exon.genome_start - 1, -1):
                    current_pos = transcript_pos + offset
                    self.genome_to_transcript[genome_pos] = current_pos
                    self.transcript_to_genome[current_pos] = genome_pos
                    offset += 1
                transcript_pos += offset

    def _check_and_extend_for_uorf(self):
        """
        Check if uORF is outside of transcript boundaries (in genomic coordinates)
        and extend the transcript if necessary.
        """
        import logging
        
        print(f"\nDEBUG: Checking if transcript {self.transcript_id} needs extension")
        print(f"uORF genomic coords: {self.uorf_start_genomic}-{self.uorf_end_genomic}")
        
        if self.uorf_start_genomic is None or self.uorf_end_genomic is None:
            print("No uORF genomic coordinates, skipping extension")
            return
        
        # Check if genome_to_transcript is empty (no exons mapped)
        if not self.genome_to_transcript:
            print(f"WARNING: Empty coordinate map for transcript {self.transcript_id}, cannot extend")
            logging.warning(f"Empty coordinate map for transcript {self.transcript_id}, cannot extend")
            return
        
        # Check if uORF is outside transcript boundaries
        if self.strand == '+':
            # For positive strand
            min_genome_pos = min(self.genome_to_transcript.keys())
            max_genome_pos = max(self.genome_to_transcript.keys())
            
            print(f"Transcript genomic range: {min_genome_pos}-{max_genome_pos} (positive strand)")
            
            # Check if uORF start is before transcript start
            if self.uorf_start_genomic < min_genome_pos:
                print(f"uORF starts before transcript {self.transcript_id} (pos: {self.uorf_start_genomic} < {min_genome_pos})")
                logging.info(f"uORF starts before transcript {self.transcript_id} (pos: {self.uorf_start_genomic} < {min_genome_pos})")
                self._extend_transcript_upstream(min_genome_pos, self.uorf_start_genomic)
                
            # Check if uORF end is after transcript end
            if self.uorf_end_genomic > max_genome_pos:
                print(f"uORF ends after transcript {self.transcript_id} (pos: {self.uorf_end_genomic} > {max_genome_pos})")
                logging.info(f"uORF ends after transcript {self.transcript_id} (pos: {self.uorf_end_genomic} > {max_genome_pos})")
                self._extend_transcript_downstream(max_genome_pos, self.uorf_end_genomic)
                
        else:  # Negative strand
            # For negative strand, the coordinate logic is reversed in genomic space
            # but the extension logic remains the same in terms of transcript positions
            min_genome_pos = min(self.genome_to_transcript.keys())
            max_genome_pos = max(self.genome_to_transcript.keys())
            
            print(f"Transcript genomic range: {min_genome_pos}-{max_genome_pos} (negative strand)")
            
            # On the negative strand, 5' end (transcript start) corresponds to max_genome_pos
            # 3' end (transcript end) corresponds to min_genome_pos
            
            # For negative strand, uORF genomic coordinates can be in reverse order
            # relative to the direction of transcription, so we check both cases
            
            # Check if uORF extends beyond the 5' end of transcript (greater than max_genome_pos)
            if self.uorf_start_genomic > max_genome_pos or self.uorf_end_genomic > max_genome_pos:
                print(f"uORF extends beyond 5' end of transcript {self.transcript_id} on negative strand")
                extend_pos = max(self.uorf_start_genomic, self.uorf_end_genomic)
                logging.info(f"Extending transcript {self.transcript_id} at 5' end (negative strand): {max_genome_pos} -> {extend_pos}")
                self._extend_transcript_downstream(max_genome_pos, extend_pos)
                
            # Check if uORF extends beyond the 3' end of transcript (less than min_genome_pos)
            if self.uorf_start_genomic < min_genome_pos or self.uorf_end_genomic < min_genome_pos:
                print(f"uORF extends beyond 3' end of transcript {self.transcript_id} on negative strand")
                extend_pos = min(self.uorf_start_genomic, self.uorf_end_genomic)
                logging.info(f"Extending transcript {self.transcript_id} at 3' end (negative strand): {min_genome_pos} -> {extend_pos}")
                self._extend_transcript_upstream(min_genome_pos, extend_pos)
                
        if self.was_extended:
            print(f"After extension: Transcript now spans {min(self.genome_to_transcript.keys())}-{max(self.genome_to_transcript.keys())}")
        else:
            print("No extension was needed")

    def _extend_transcript_upstream(self, current_min_pos, new_min_pos):
        """
        Extend transcript upstream (towards 5' end depending on strand).
        
        Args:
            current_min_pos: Current minimum genomic position in transcript
            new_min_pos: New minimum genomic position to extend to
        """
        import logging
        
        print(f"Extending transcript {self.transcript_id} upstream: {current_min_pos} -> {new_min_pos}")
        
        # Ensure proper ordering
        if new_min_pos > current_min_pos:
            print(f"WARNING: Cannot extend upstream, new position {new_min_pos} > current position {current_min_pos}")
            return
            
        # Extension length
        extension_length = current_min_pos - new_min_pos
        print(f"Extension length: {extension_length}")
        
        # Create temporary maps for the extension
        new_genome_to_transcript = {}
        new_transcript_to_genome = {}
        
        # Shift all existing transcript positions by extension length
        for genome_pos, transcript_pos in self.genome_to_transcript.items():
            new_genome_to_transcript[genome_pos] = transcript_pos + extension_length
            
        for transcript_pos, genome_pos in self.transcript_to_genome.items():
            new_transcript_to_genome[transcript_pos + extension_length] = genome_pos
        
        print(f"Shifted {len(self.genome_to_transcript)} existing positions by {extension_length}")
        
        # Add new mappings for extended region
        for i in range(extension_length):
            genome_pos = new_min_pos + i
            transcript_pos = i + 1  # Start from 1 (1-based)
            new_genome_to_transcript[genome_pos] = transcript_pos
            new_transcript_to_genome[transcript_pos] = genome_pos
        
        print(f"Added {extension_length} new positions to map")
        
        # Update maps
        self.genome_to_transcript = new_genome_to_transcript
        self.transcript_to_genome = new_transcript_to_genome
        
        self.was_extended = True
        print(f"Transcript extension complete. Min position: {min(self.genome_to_transcript.keys())}")
        logging.info(f"Extended transcript {self.transcript_id} upstream by {extension_length} bp")

    def _extend_transcript_downstream(self, current_max_pos, new_max_pos):
        """
        Extend transcript downstream (towards 3' end depending on strand).
        
        Args:
            current_max_pos: Current maximum genomic position in transcript
            new_max_pos: New maximum genomic position to extend to
        """
        import logging
        
        print(f"Extending transcript {self.transcript_id} downstream: {current_max_pos} -> {new_max_pos}")
        
        # Ensure proper ordering
        if new_max_pos < current_max_pos:
            print(f"WARNING: Cannot extend downstream, new position {new_max_pos} < current position {current_max_pos}")
            return
            
        # Extension length
        extension_length = new_max_pos - current_max_pos
        print(f"Extension length: {extension_length}")
        
        # Get the maximum transcript position
        max_transcript_pos = max(self.transcript_to_genome.keys()) if self.transcript_to_genome else 0
        print(f"Current max transcript position: {max_transcript_pos}")
        
        # Add new mappings for extended region
        count = 0
        for i in range(extension_length):
            genome_pos = current_max_pos + i + 1
            transcript_pos = max_transcript_pos + i + 1
            self.genome_to_transcript[genome_pos] = transcript_pos
            self.transcript_to_genome[transcript_pos] = genome_pos
            count += 1
        
        print(f"Added {count} new positions to map")
        
        self.was_extended = True
        print(f"Transcript extension complete. Max position: {max(self.genome_to_transcript.keys())}")
        logging.info(f"Extended transcript {self.transcript_id} downstream by {extension_length} bp")

    def get_coordinates(self, genomic_pos: int) -> Optional[Coordinates]:
        """Get both genomic and transcript coordinates for a position."""
        transcript_pos = self.genome_to_transcript.get(genomic_pos)
        if transcript_pos is not None:
            return Coordinates(genomic_pos, transcript_pos)
        return None

    def get_genomic_position(self, transcript_pos: int) -> Optional[int]:
        """Convert transcript position to genomic position."""
        return self.transcript_to_genome.get(transcript_pos)

    def get_transcript_position(self, genomic_pos: int) -> Optional[int]:
        """Convert genomic position to transcript position."""
        return self.genome_to_transcript.get(genomic_pos)

    @property
    def mainorf_coords(self) -> Optional[Coordinates]:
        """Get main ORF coordinates."""
        if self.mainorf_start_genomic and self.mainorf_start is not None:
            return Coordinates(self.mainorf_start_genomic, self.mainorf_start)
        return None

    @property
    def uorf_coords(self) -> Optional[Coordinates]:
        """Get uORF coordinates."""
        if self.uorf_start_genomic and self.uorf_start is not None:
            return Coordinates(self.uorf_start_genomic, self.uorf_start)
        return None