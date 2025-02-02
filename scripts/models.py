# models.py

from dataclasses import dataclass, field
from typing import List, Dict, Optional

@dataclass
class Exon:
    """
    Representation of an exon with genomic and transcript coordinates.
    
    Attributes:
        start: Starting position in transcript coordinates
        length: Length of the exon
        genome_start: Starting position in genomic coordinates (1-based)
        genome_end: Ending position in genomic coordinates (1-based)
    """
    start: int
    length: int
    genome_start: int
    genome_end: int

    def __post_init__(self):
        """Validate exon coordinates after initialization."""
        if self.genome_start > self.genome_end:
            raise ValueError(f"Invalid exon coordinates: start ({self.genome_start}) > end ({self.genome_end})")
        if self.length != (self.genome_end - self.genome_start + 1):
            raise ValueError(f"Exon length mismatch: {self.length} != {self.genome_end - self.genome_start + 1}")

@dataclass
class Transcript:
    """
    Representation of a transcript with its exons and coordinate mappings.
    
    Attributes:
        transcript_id: Unique identifier for the transcript
        chromosome: Chromosome name
        strand: Strand ('+' or '-')
        exons: List of Exon objects
        genome_to_transcript: Mapping from genomic to transcript coordinates
        transcript_to_genome: Mapping from transcript to genomic coordinates
    """
    transcript_id: str
    chromosome: str
    strand: str
    exons: List[Exon]
    genome_to_transcript: Dict[int, int] = field(default_factory=dict)
    transcript_to_genome: Dict[int, int] = field(default_factory=dict)
    
    def __post_init__(self):
        """Initialize transcript after creation."""
        self._validate_strand()
        self._sort_exons()
        self._build_coordinate_maps()
    
    def _validate_strand(self):
        """Validate strand information."""
        if self.strand not in ['+', '-']:
            raise ValueError(f"Invalid strand: {self.strand}. Must be '+' or '-'")
    
    def _sort_exons(self):
        """Sort exons by genomic start position."""
        self.exons = sorted(self.exons, key=lambda x: x.genome_start)
        
        # Check for overlapping exons
        for i in range(len(self.exons) - 1):
            if self.exons[i].genome_end >= self.exons[i + 1].genome_start:
                raise ValueError(f"Overlapping exons detected in transcript {self.transcript_id}")
    
    def _build_coordinate_maps(self):
        """Build mappings between genomic and transcript coordinates."""
        self.genome_to_transcript.clear()
        self.transcript_to_genome.clear()
        
        if self.strand == '+':
            current_pos = 1  # 1-based coordinates
            for exon in self.exons:
                for genome_pos in range(exon.genome_start, exon.genome_end + 1):
                    self.genome_to_transcript[genome_pos] = current_pos
                    self.transcript_to_genome[current_pos] = genome_pos
                    current_pos += 1
        else:
            current_pos = 1
            for exon in reversed(self.exons):
                for genome_pos in range(exon.genome_end, exon.genome_start - 1, -1):
                    self.genome_to_transcript[genome_pos] = current_pos
                    self.transcript_to_genome[current_pos] = genome_pos
                    current_pos += 1
    
    def get_exon_by_position(self, genome_pos: int) -> Optional[Exon]:
        """
        Find exon containing the given genomic position.
        
        Args:
            genome_pos: Genomic position
            
        Returns:
            Exon object if position is in an exon, None otherwise
        """
        for exon in self.exons:
            if exon.genome_start <= genome_pos <= exon.genome_end:
                return exon
        return None
    
    def get_total_transcript_length(self) -> int:
        """Calculate total transcript length."""
        return sum(exon.length for exon in self.exons)
    
    def is_position_in_transcript(self, genome_pos: int) -> bool:
        """Check if a genomic position falls within the transcript."""
        return genome_pos in self.genome_to_transcript