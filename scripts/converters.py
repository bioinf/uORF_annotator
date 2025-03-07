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
        
        print("\n=== DEBUG: BED Entries ===")
        for idx, bed_row in bed_df.iterrows():
            print(f"BED Row {idx}: {bed_row[0]}, {bed_row[1]}, {bed_row[2]}, {bed_row[3]}, strand: {bed_row[5]}")
        print("=== END BED Entries ===\n")
        
        # Create a dictionary to track which transcripts have been processed
        processed_transcripts = {}
        
        for _, bed_row in bed_df.iterrows():
            transcript_id = bed_row[3].split('|')[0].split('.')[0]
            uorf_id = f"{transcript_id}_uorf_{len(processed_transcripts.get(transcript_id, []))+1}"
            
            print(f"\nProcessing transcript: {transcript_id}, uORF ID: {uorf_id}")
            print(f"uORF genomic coords (from BED): {int(bed_row[1]) + 1}-{int(bed_row[2])}")
            
            if transcript_id in gtf_data['exons']:
                print(f"Found exons for {transcript_id}: {len(gtf_data['exons'][transcript_id])} exons")
                mainorf_coords = self._get_mainorf_coords(transcript_id, gtf_data['cds'])
                print(f"MainORF genomic coords: {mainorf_coords[0]}-{mainorf_coords[1]}")

                # Check if we already created a transcript object for this transcript_id
                if transcript_id in self.transcripts:
                    print(f"Transcript {transcript_id} already exists, adding another uORF")
                    
                    # Get the existing transcript
                    existing_transcript = self.transcripts[transcript_id]
                    
                    # Create a new transcript with the same properties but different uORF
                    new_transcript = Transcript(
                        transcript_id=uorf_id,  # Use unique ID for the new transcript-uORF pair
                        chromosome=existing_transcript.chromosome,
                        strand=existing_transcript.strand,
                        exons=existing_transcript.exons.copy(),
                        mainorf_start=existing_transcript.mainorf_start_genomic,
                        mainorf_end=existing_transcript.mainorf_end_genomic,
                        uorf_start=int(bed_row[1]) + 1,
                        uorf_end=int(bed_row[2])
                    )
                    
                    # Store the new transcript-uORF pair
                    self.transcripts[uorf_id] = new_transcript
                    
                    # Track that we've processed this transcript
                    if transcript_id not in processed_transcripts:
                        processed_transcripts[transcript_id] = []
                    processed_transcripts[transcript_id].append(uorf_id)
                    
                    # Debugging output
                    t = new_transcript
                    print(f"Created additional uORF transcript {uorf_id} for {transcript_id}:")
                else:
                    # Create the first transcript-uORF pair for this transcript_id
                    self.transcripts[transcript_id] = Transcript(
                        transcript_id=transcript_id,
                        chromosome=bed_row[0],
                        strand=bed_row[5],
                        exons=gtf_data['exons'][transcript_id],
                        mainorf_start=mainorf_coords[0],
                        mainorf_end=mainorf_coords[1],
                        uorf_start=int(bed_row[1]) + 1,
                        uorf_end=int(bed_row[2])
                    )
                    
                    # Track that we've processed this transcript
                    processed_transcripts[transcript_id] = [transcript_id]
                    
                    # Debugging output
                    t = self.transcripts[transcript_id]
                    print(f"Created transcript {transcript_id}:")
                
                # Common debugging for both cases
                print(f"  Strand: {t.strand}")
                print(f"  uORF genomic: {t.uorf_start_genomic}-{t.uorf_end_genomic}")
                print(f"  uORF transcript: {t.uorf_start}-{t.uorf_end}")
                print(f"  MainORF genomic: {t.mainorf_start_genomic}-{t.mainorf_end_genomic}")
                print(f"  MainORF transcript: {t.mainorf_start}-{t.mainorf_end}")
                print(f"  Was extended: {getattr(t, 'was_extended', False)}")
                
                # Checks for genomic range
                min_g = min(t.genome_to_transcript.keys()) if t.genome_to_transcript else None
                max_g = max(t.genome_to_transcript.keys()) if t.genome_to_transcript else None
                print(f"  Genomic range: {min_g}-{max_g}")
            else:
                print(f"WARNING: No exons found for transcript {transcript_id}")

        # Print summary of transcripts with multiple uORFs
        print("\n=== Multiple uORF Summary ===")
        for transcript_id, uorf_ids in processed_transcripts.items():
            if len(uorf_ids) > 1:
                print(f"Transcript {transcript_id} has {len(uorf_ids)} uORFs: {', '.join(uorf_ids)}")
        print("===========================\n")

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