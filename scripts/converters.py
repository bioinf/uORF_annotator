"""
Coordinate conversion utilities for genomic and transcript coordinates.
"""

import logging
import gzip
from typing import Dict, Tuple
import pandas as pd

from models import Exon, Transcript
from parsers import GTFParser


class CoordinateConverter:
    """
    Handles conversion between genomic and transcript coordinates.
    Manages transcript information and coordinate mappings.
    """
    
    def __init__(self, bed_file: str, gtf_file: str):
        """
        Initialize the converter with input files.
        
        Args:
            bed_file: Path to BED format file
            gtf_file: Path to GTF format file
        """
        self.bed_file = bed_file
        self.gtf_file = gtf_file
        self.transcripts: Dict[str, Transcript] = {}
        self._parse_files()

    def _parse_files(self) -> None:
        """Parse both BED and GTF files to build transcript information."""
        gtf_data = self._parse_gtf()
        self._parse_bed(gtf_data)

    def _parse_gtf(self) -> Dict:
        """
        Parse GTF file to extract CDS coordinates.
        
        Returns:
            Dictionary containing CDS and exon information
        """
        cds_by_transcript = {}
        exons_by_transcript = {}
        
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parsed = GTFParser.parse_line(line)
                if not parsed:
                    continue

                transcript_id = parsed['transcript_id']
                if not transcript_id:
                    continue

                self._process_gtf_feature(parsed, transcript_id, 
                                       exons_by_transcript, cds_by_transcript)

        return {'cds': cds_by_transcript, 'exons': exons_by_transcript}

    def _process_gtf_feature(self, parsed: dict, transcript_id: str,
                           exons_by_transcript: Dict, cds_by_transcript: Dict) -> None:
        """
        Process individual features from GTF file.
        
        Args:
            parsed: Parsed GTF line data
            transcript_id: ID of the transcript
            exons_by_transcript: Dictionary to store exon information
            cds_by_transcript: Dictionary to store CDS information
        """
        if parsed['feature'] == 'exon':
            if transcript_id not in exons_by_transcript:
                exons_by_transcript[transcript_id] = []
            exon = Exon(
                start=len(exons_by_transcript[transcript_id]),
                length=parsed['end'] - parsed['start'] + 1,
                genome_start=parsed['start'],
                genome_end=parsed['end']
            )
            exons_by_transcript[transcript_id].append(exon)
        
        elif parsed['feature'] == 'CDS':
            if transcript_id not in cds_by_transcript:
                cds_by_transcript[transcript_id] = {'starts': [], 'ends': []}
            cds_by_transcript[transcript_id]['starts'].append(parsed['start'])
            cds_by_transcript[transcript_id]['ends'].append(parsed['end'])

    def _parse_bed(self, gtf_data: Dict) -> None:
        """
        Parse BED file and create transcript objects.
        
        Args:
            gtf_data: Dictionary containing CDS and exon information from GTF
        """
        bed_df = pd.read_csv(self.bed_file, sep='\t', header=None)
        
        for _, bed_row in bed_df.iterrows():
            transcript_id = bed_row[3].split('|')[0].split('.')[0]
            
            if transcript_id in gtf_data['exons']:
                mainorf_coords = self._get_mainorf_coords(transcript_id, gtf_data['cds'])

                self.transcripts[transcript_id] = Transcript(
                    transcript_id=transcript_id,
                    chromosome=bed_row[0],
                    strand=bed_row[5],
                    exons=gtf_data['exons'][transcript_id],
                    mainorf_start=mainorf_coords[0],
                    mainorf_end=mainorf_coords[1]
                )

    def _get_mainorf_coords(self, transcript_id: str, cds_data: Dict) -> Tuple[int, int]:
        """
        Get main ORF coordinates for a transcript.
        
        Args:
            transcript_id: ID of the transcript
            cds_data: Dictionary containing CDS information
            
        Returns:
            Tuple of (start, end) coordinates for main ORF
        """
        if transcript_id in cds_data:
            return (
                min(cds_data[transcript_id]['starts']),
                max(cds_data[transcript_id]['ends'])
            )
        return None, None

    def genome_to_transcript_pos(self, transcript_id: str, genome_pos: int) -> str:
        """
        Convert genomic position to transcript position.
        
        Args:
            transcript_id: ID of the transcript
            genome_pos: Position in genomic coordinates
            
        Returns:
            Position in transcript coordinates or "NA" if conversion fails
        """
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            logging.warning(f"Transcript {transcript_id} not found")
            return "NA"

        return transcript.genome_to_transcript.get(genome_pos, "NA")