"""
Models for representing genomic and transcript data structures.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional


@dataclass
class Exon:
    """
    Represents an exon with both genomic and transcript coordinates.
    
    Attributes:
        start: Starting position in the transcript
        length: Length of the exon
        genome_start: Starting position in the genome
        genome_end: Ending position in the genome
    """
    start: int
    length: int
    genome_start: int
    genome_end: int


class Transcript:
    """
    Represents a transcript with methods for coordinate conversion.
    
    Handles both forward and reverse strand transcripts, maintaining
    mappings between genomic and transcript coordinates.
    """
    
    def __init__(self, transcript_id: str, chromosome: str, strand: str, 
                 exons: List[Exon], mainorf_start: Optional[int] = None, 
                 mainorf_end: Optional[int] = None):
        """
        Initialize a transcript with its basic properties and exon structure.
        
        Args:
            transcript_id: Unique identifier for the transcript
            chromosome: Chromosome name/number
            strand: DNA strand ('+' or '-')
            exons: List of Exon objects
            mainorf_start: Start position of main ORF (optional)
            mainorf_end: End position of main ORF (optional)
        """
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = sorted(exons, key=lambda x: x.genome_start)
        self.mainorf_start = mainorf_start
        self.mainorf_end = mainorf_end
        self.genome_to_transcript = {}
        self.transcript_to_genome = {}
        self._build_coordinate_maps()

    def _build_coordinate_maps(self):
        current_transcript_pos = 0
        if self.strand == '+':
            for exon in self.exons:
                # print(f"Processing forward strand exon: {exon.genome_start}-{exon.genome_end}")
                for genome_pos in range(exon.genome_start, exon.genome_end + 1):
                    self.genome_to_transcript[genome_pos] = current_transcript_pos
                    self.transcript_to_genome[current_transcript_pos] = genome_pos
                    current_transcript_pos += 1
        else:
            print(f"Building coordinates for reverse strand transcript {self.transcript_id}")
            for exon in reversed(self.exons):
                # print(f"Processing reverse strand exon: {exon.genome_start}-{exon.genome_end}")
                for genome_pos in range(exon.genome_end, exon.genome_start - 1, -1):
                    # print(f"Mapping genome pos {genome_pos} to transcript pos {current_transcript_pos}")
                    self.genome_to_transcript[genome_pos] = current_transcript_pos
                    self.transcript_to_genome[current_transcript_pos] = genome_pos
                    current_transcript_pos += 1

    def _add_position_mapping(self, genome_pos: int, transcript_pos: int) -> None:
        """
        Add a bidirectional mapping between genomic and transcript positions.
        
        Args:
            genome_pos: Position in genomic coordinates
            transcript_pos: Position in transcript coordinates
        """
        self.genome_to_transcript[genome_pos] = transcript_pos
        self.transcript_to_genome[transcript_pos] = genome_pos