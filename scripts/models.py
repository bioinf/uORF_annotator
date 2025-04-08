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
                 uorf_end: Optional[int] = None,
                 start_codon_type: Optional[str] = None,  # New parameter to store ATG vs non-ATG status
                 debug_mode: bool = False):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = sorted(exons, key=lambda x: x.genome_start)
        self.debug_mode = debug_mode
        
        # Store genomic coordinates
        self.mainorf_start_genomic = mainorf_start
        self.mainorf_end_genomic = mainorf_end
        self.uorf_start_genomic = uorf_start
        self.uorf_end_genomic = uorf_end
        
        # Store start codon type (ATG, NON-ATG, or None if unknown)
        self.start_codon_type = start_codon_type
        
        # Initialize coordinate maps
        self.genome_to_transcript = {}
        self.transcript_to_genome = {}
        
        # Flag to indicate if the transcript was extended to include uORF
        self.was_extended = False
        
        # Flag for overlap between uORF and mainCDS (will be set by CoordinateConverter)
        self.overlaps_maincds = False
        
        # Store original coordinates before any manipulation
        self._store_original_coordinates()
        
        # Build coordinate map for the transcript, handling extension if needed
        self._rebuild_coordinates_for_transcript()
        
        # Ensure coordinates are consistent with strand orientation
        self._normalize_coordinates_by_strand()

    def _store_original_coordinates(self):
        """Store original coordinates before any manipulation."""
        # Store original genomic coordinates
        self.original_mainorf_start_genomic = self.mainorf_start_genomic
        self.original_mainorf_end_genomic = self.mainorf_end_genomic
        self.original_uorf_start_genomic = self.uorf_start_genomic
        self.original_uorf_end_genomic = self.uorf_end_genomic

    def _convert_genomic_to_transcript_coords(self):
        """Convert all genomic coordinates to transcript coordinates."""
        self.mainorf_start = self._convert_to_transcript_coords(self.mainorf_start_genomic)
        self.mainorf_end = self._convert_to_transcript_coords(self.mainorf_end_genomic)
        self.uorf_start = self._convert_to_transcript_coords(self.uorf_start_genomic)
        self.uorf_end = self._convert_to_transcript_coords(self.uorf_end_genomic)
        
        # Store original transcript coordinates
        self.original_mainorf_start = self.mainorf_start
        self.original_mainorf_end = self.mainorf_end
        self.original_uorf_start = self.uorf_start
        self.original_uorf_end = self.uorf_end

    def _normalize_coordinates_by_strand(self):
        """
        Ensure coordinates are properly ordered based on strand.
        For positive strand: start < end
        For negative strand, we maintain biological ordering in transcript coordinates,
        which is already handled by the coordinate conversion process.
        """
        # Check if any coordinates are missing
        if self.uorf_start is None or self.uorf_end is None:
            logging.warning(f"Missing uORF transcript coordinates for {self.transcript_id}")
            return
            
        if self.mainorf_start is None or self.mainorf_end is None:
            logging.warning(f"Missing mainORF transcript coordinates for {self.transcript_id}")
            return
            
        # For positive strand, ensure start < end
        if self.strand == '+':
            # Fix uORF coordinates if needed
            if self.uorf_start > self.uorf_end:
                logging.warning(f"Fixing uORF coordinates for positive strand: {self.uorf_start} > {self.uorf_end}")
                self.uorf_start, self.uorf_end = self.uorf_end, self.uorf_start
                
            # Fix mainORF coordinates if needed
            if self.mainorf_start > self.mainorf_end:
                logging.warning(f"Fixing mainORF coordinates for positive strand: {self.mainorf_start} > {self.mainorf_end}")
                self.mainorf_start, self.mainorf_end = self.mainorf_end, self.mainorf_start
        
        # For negative strand, we might need the opposite ordering in some cases
        # but this depends on the specific use case and is handled elsewhere if needed

    def _convert_to_transcript_coords(self, genomic_pos: int) -> Optional[int]:
        """
        Convert genomic position to transcript coordinates.
        
        Args:
            genomic_pos: Genomic position to convert
            
        Returns:
            Corresponding transcript position, or None if conversion fails
        """
        if genomic_pos is None:
            return None
            
        # Get transcript position from mapping
        transcript_pos = self.genome_to_transcript.get(genomic_pos)
        
        # Log debug info for extended transcripts
        if transcript_pos is None and self.debug_mode:
            logging.debug(f"No transcript position found for genomic position {genomic_pos} in transcript {self.transcript_id}")
            
            # For debugging purposes, show the range of valid genomic positions
            if self.genome_to_transcript:
                min_pos = min(self.genome_to_transcript.keys())
                max_pos = max(self.genome_to_transcript.keys())
                logging.debug(f"Valid genomic positions range: {min_pos}-{max_pos}")
        
        return transcript_pos

    def _rebuild_coordinates_for_transcript(self):
        """
        Rebuild the coordinate map for the entire transcript, ensuring
        consistent 5' to 3' ordering in transcript coordinates.
        This method replaces _build_coordinate_maps and handles extension.
        """
        # Start fresh
        self.genome_to_transcript = {}
        self.transcript_to_genome = {}
        
        # Get all exons - we'll work with these directly
        if not self.exons:
            logging.warning(f"No exons defined for transcript {self.transcript_id}")
            return
        
        # Sort exons by genomic position
        sorted_exons = sorted(self.exons, key=lambda x: x.genome_start)
        
        # Check if uORF extends beyond transcript boundaries
        needs_extension = False
        min_transcript_pos = sorted_exons[0].genome_start
        max_transcript_pos = sorted_exons[-1].genome_end
        
        # Define the full genomic range to include
        start_pos = min_transcript_pos
        end_pos = max_transcript_pos
        
        # Check if we need to extend to include uORF
        if self.uorf_start_genomic and self.uorf_end_genomic:
            uorf_min = min(self.uorf_start_genomic, self.uorf_end_genomic)
            uorf_max = max(self.uorf_start_genomic, self.uorf_end_genomic)
            
            if uorf_min < start_pos:
                start_pos = uorf_min
                needs_extension = True
                logging.info(f"Extending transcript {self.transcript_id} to include uORF start")
                
            if uorf_max > end_pos:
                end_pos = uorf_max
                needs_extension = True
                logging.info(f"Extending transcript {self.transcript_id} to include uORF end")
        
        # For negative strand, we need to reverse the exon order for transcript coordinates
        if self.strand == '-':
            # For negative strand, transcript positions start at the 5' end
            # This means we count from the highest genomic coordinate
            exon_positions = []
            for exon in sorted_exons:
                # Include only genomic ranges within exons (not introns)
                exon_positions.extend(range(exon.genome_start, exon.genome_end + 1))
            
            # If we need to extend to include uORF
            if needs_extension:
                # Add extension at 5' end (higher genomic coordinates)
                if end_pos > max_transcript_pos:
                    extension_range = range(max_transcript_pos + 1, end_pos + 1)
                    exon_positions.extend(extension_range)
                    
                # Add extension at 3' end (lower genomic coordinates)
                if start_pos < min_transcript_pos:
                    extension_range = range(start_pos, min_transcript_pos)
                    exon_positions.extend(extension_range)
            
            # Now sort all positions in descending order (for negative strand)
            # This ensures transcript positions increase from 5' to 3'
            exon_positions.sort(reverse=True)
            
            # Assign transcript positions
            transcript_pos = 1
            for genome_pos in exon_positions:
                self.genome_to_transcript[genome_pos] = transcript_pos
                self.transcript_to_genome[transcript_pos] = genome_pos
                transcript_pos += 1
                
        else:  # Positive strand
            # For positive strand, transcript positions increase with genomic positions
            exon_positions = []
            for exon in sorted_exons:
                # Include only genomic ranges within exons (not introns)
                exon_positions.extend(range(exon.genome_start, exon.genome_end + 1))
            
            # If we need to extend to include uORF
            if needs_extension:
                # Add extension at 5' end (lower genomic coordinates)
                if start_pos < min_transcript_pos:
                    extension_range = range(start_pos, min_transcript_pos)
                    exon_positions.extend(extension_range)
                    
                # Add extension at 3' end (higher genomic coordinates)
                if end_pos > max_transcript_pos:
                    extension_range = range(max_transcript_pos + 1, end_pos + 1)
                    exon_positions.extend(extension_range)
            
            # Sort all positions
            exon_positions.sort()
            
            # Assign transcript positions
            transcript_pos = 1
            for genome_pos in exon_positions:
                self.genome_to_transcript[genome_pos] = transcript_pos
                self.transcript_to_genome[transcript_pos] = genome_pos
                transcript_pos += 1
        
        self.was_extended = needs_extension
        logging.info(f"Rebuilt coordinate map for transcript {self.transcript_id} ({self.strand} strand), "
                    f"length: {len(self.genome_to_transcript)} bp, extended: {self.was_extended}")

        # After building coordinates, convert genomic to transcript coordinates
        self._convert_genomic_to_transcript_coords()

    def get_coordinates(self, genomic_pos: int) -> Optional[Coordinates]:
        """Get both genomic and transcript coordinates for a position with detailed logging."""
        transcript_pos = self.genome_to_transcript.get(genomic_pos)
        if transcript_pos is not None:
            if self.debug_mode:
                logging.debug(f"Genomic pos {genomic_pos} -> Transcript pos {transcript_pos} "
                        f"for {self.transcript_id} (strand {self.strand}, extended: {self.was_extended})")
            return Coordinates(genomic_pos, transcript_pos)
        
        logging.warning(f"Genomic position {genomic_pos} not found in coordinate map for "
                      f"transcript {self.transcript_id} (ext: {self.was_extended})")
                      
        # For debugging, show the range of valid positions
        if self.genome_to_transcript and self.debug_mode:
            min_pos = min(self.genome_to_transcript.keys())
            max_pos = max(self.genome_to_transcript.keys())
            logging.warning(f"Valid genomic positions range: {min_pos}-{max_pos}")
        
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