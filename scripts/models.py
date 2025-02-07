"""
Models for representing genomic and transcript data structures.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, NamedTuple


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
        
        # Build coordinate maps
        self._build_coordinate_maps()
        
        # Convert and store transcript coordinates
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