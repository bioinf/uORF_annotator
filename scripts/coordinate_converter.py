# coordinate_converter.py

from tqdm import tqdm
import logging
import gzip
from typing import Dict, Optional, Union, List, Tuple
import pandas as pd
from models import Exon, Transcript

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CoordinateConverter:
    """Convert between genomic and transcript coordinates."""
    
    def __init__(self, bed_file: str, gtf_file: str):
        """
        Initialize converter with input files.
        
        Args:
            bed_file: Path to BED file with uORF coordinates
            gtf_file: Path to GTF file with genomic features
        """
        self.bed_file = bed_file
        self.gtf_file = gtf_file
        self.transcripts: Dict[str, Transcript] = {}
        self.uorf_coords: Dict[str, Dict[str, int]] = {}  # transcript_id -> {'start': x, 'end': y}
        self.cds_coords: Dict[str, Dict[str, int]] = {}   # transcript_id -> {'start': x, 'end': y}
        self._parse_files()
        
    def _parse_bed_line(self, row: pd.Series) -> Optional[dict]:
        """Parse a single BED line for uORF coordinates."""
        try:
            transcript_id = row[3].split('|')[0].split('.')[0]  # Get NM_XXXXX
            return {
                'transcript_id': transcript_id,
                'chromosome': row[0],
                'start': int(row[1]),
                'end': int(row[2]),
                'strand': row[5]
            }
        except Exception as e:
            logger.error(f"Error parsing BED line: {e}")
            return None

    def _parse_gtf_line(self, line: str) -> Optional[dict]:
        """
        Parse a single GTF line and return relevant information.
        
        Args:
            line: Single line from GTF file
            
        Returns:
            Dictionary with parsed information or None if parsing fails
        """
        try:
            fields = line.strip().split('\t')
            if len(fields) < 9:
                return None

            # Parse attributes
            attributes = {}
            for attr in fields[8].rstrip(';').split('; '):
                if ' ' in attr:
                    key, value = attr.strip().split(' ', 1)
                    attributes[key] = value.strip('"')

            transcript_id = attributes.get('transcript_id', '')
            if not transcript_id:
                return None

            return {
                'chromosome': fields[0],
                'feature': fields[2],
                'start': int(fields[3]),  # GTF is 1-based
                'end': int(fields[4]),
                'strand': fields[6],
                'transcript_id': transcript_id.split('.')[0]  # Remove version
            }
        except Exception as e:
            logger.error(f"Error parsing GTF line: {e}")
            return None

    def _parse_files(self):
        """Parse BED and GTF files to build transcript information."""
        logger.info("Starting to parse input files...")

        # Read BED file
        logger.info("Reading BED file...")
        bed_df = pd.read_csv(self.bed_file, sep='\t', header=None)

        # Process BED file for uORF coordinates
        for _, row in bed_df.iterrows():
            bed_info = self._parse_bed_line(row)
            if bed_info:
                self.uorf_coords[bed_info['transcript_id']] = {
                    'start': bed_info['start'],
                    'end': bed_info['end']
                }

        # Initialize dictionary to store exons by transcript
        exons_by_transcript: Dict[str, List[Exon]] = {}

        # Parse GTF file
        logger.info("Reading GTF file...")
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        line_count = 0

        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                line_count += 1
                if line_count % 10000 == 0:
                    logger.info(f"Processed {line_count} lines from GTF")

                if line.startswith('#'):
                    continue

                parsed = self._parse_gtf_line(line)
                if not parsed:
                    continue

                transcript_id = parsed['transcript_id']

                # Process CDS feature
                if parsed['feature'] == 'CDS':
                    if transcript_id not in self.cds_coords:
                        self.cds_coords[transcript_id] = {
                            'start': parsed['start'],
                            'end': parsed['start']  # Initialize with start
                        }
                    # Update CDS boundaries
                    self.cds_coords[transcript_id]['start'] = min(
                        self.cds_coords[transcript_id]['start'], 
                        parsed['start']
                    )
                    self.cds_coords[transcript_id]['end'] = max(
                        self.cds_coords[transcript_id]['end'], 
                        parsed['end']
                    )

                # Process exon feature
                elif parsed['feature'] == 'exon':
                    if transcript_id not in exons_by_transcript:
                        exons_by_transcript[transcript_id] = []

                    exon = Exon(
                        start=len(exons_by_transcript[transcript_id]),
                        length=parsed['end'] - parsed['start'] + 1,
                        genome_start=parsed['start'],
                        genome_end=parsed['end']
                    )
                    exons_by_transcript[transcript_id].append(exon)

        # Create Transcript objects
        logger.info("Creating Transcript objects...")
        for transcript_id in tqdm(exons_by_transcript.keys(), desc="Creating transcripts"):
            if transcript_id in self.uorf_coords:  # Only create transcripts that have uORFs
                try:
                    # Find BED entry for this transcript
                    bed_info = None
                    for _, row in bed_df.iterrows():
                        parsed_bed = self._parse_bed_line(row)
                        if parsed_bed and parsed_bed['transcript_id'] == transcript_id:
                            bed_info = parsed_bed
                            break
                    
                    if bed_info and transcript_id in exons_by_transcript:
                        self.transcripts[transcript_id] = Transcript(
                            transcript_id=transcript_id,
                            chromosome=bed_info['chromosome'],
                            strand=bed_info['strand'],
                            exons=exons_by_transcript[transcript_id]
                        )
                except Exception as e:
                    logger.warning(f"Skipping transcript {transcript_id}: {str(e)}")

        logger.info(f"Successfully created {len(self.transcripts)} transcript objects")
        logger.info(f"Found CDS coordinates for {len(self.cds_coords)} transcripts")
        logger.info(f"Found uORF coordinates for {len(self.uorf_coords)} transcripts")

    def genome_to_transcript_pos(self, transcript_id: str, genome_pos: int) -> Union[int, str]:
        """
        Convert genomic position to transcript position.
        
        Args:
            transcript_id: Transcript identifier
            genome_pos: Genomic position (1-based)
            
        Returns:
            Transcript position (1-based) or "NA" if conversion fails
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            logger.warning(f"Transcript {transcript_id} not found")
            return "NA"

        return transcript.genome_to_transcript.get(genome_pos, "NA")

    def transcript_to_genome_pos(self, transcript_id: str, transcript_pos: int) -> Union[int, str]:
        """
        Convert transcript position to genomic position.
        
        Args:
            transcript_id: Transcript identifier
            transcript_pos: Transcript position (1-based)
            
        Returns:
            Genomic position (1-based) or "NA" if conversion fails
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            logger.warning(f"Transcript {transcript_id} not found")
            return "NA"

        return transcript.transcript_to_genome.get(transcript_pos, "NA")

    def get_transcript_length(self, transcript_id: str) -> Optional[int]:
        """
        Get the total length of a transcript.
        
        Args:
            transcript_id: Transcript identifier
            
        Returns:
            Total length of transcript or None if transcript not found
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            return None
        return transcript.get_total_transcript_length()

    def get_transcript_info(self, transcript_id: str) -> Optional[dict]:
        """
        Get information about a transcript.
        
        Args:
            transcript_id: Transcript identifier
            
        Returns:
            Dictionary with transcript information or None if not found
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            return None
            
        return {
            'transcript_id': transcript.transcript_id,
            'chromosome': transcript.chromosome,
            'strand': transcript.strand,
            'length': transcript.get_total_transcript_length(),
            'num_exons': len(transcript.exons),
            'genomic_span': transcript.exons[-1].genome_end - transcript.exons[0].genome_start + 1
        }

    def check_position_in_exons(self, transcript_id: str, genome_pos: int) -> bool:
        """
        Check if a genomic position falls within any exon of the transcript.
        
        Args:
            transcript_id: Transcript identifier
            genome_pos: Genomic position to check
            
        Returns:
            True if position is in an exon, False otherwise
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            return False
            
        return transcript.is_position_in_transcript(genome_pos)

    def get_nearest_exon_boundaries(self, transcript_id: str, genome_pos: int) -> Tuple[Optional[int], Optional[int]]:
        """
        Find the nearest exon boundaries for a given genomic position.
        
        Args:
            transcript_id: Transcript identifier
            genome_pos: Genomic position
            
        Returns:
            Tuple of (previous_exon_end, next_exon_start) or (None, None) if not found
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            return None, None

        exon = transcript.get_exon_by_position(genome_pos)
        if exon:
            return exon.genome_start, exon.genome_end

        # Find nearest exon boundaries
        prev_end = None
        next_start = None

        for i, exon in enumerate(transcript.exons):
            if genome_pos < exon.genome_start:
                next_start = exon.genome_start
                if i > 0:
                    prev_end = transcript.exons[i-1].genome_end
                break
            elif genome_pos > exon.genome_end:
                prev_end = exon.genome_end
                if i < len(transcript.exons) - 1:
                    next_start = transcript.exons[i+1].genome_start

        return prev_end, next_start