"""
Parsers for handling GTF and other genomic file formats.
"""

import logging
from typing import Optional, Dict


class GTFParser:
    """Parser for GTF (Gene Transfer Format) files."""
    
    @staticmethod
    def parse_line(line: str) -> Optional[dict]:
        """
        Parse a single line from a GTF file.
        
        Args:
            line: A single line from a GTF file
            
        Returns:
            Dictionary containing parsed features or None if parsing fails
            
        The returned dictionary contains:
            - chromosome: Chromosome name
            - feature: Feature type (e.g., exon, CDS)
            - start: Feature start position
            - end: Feature end position
            - strand: DNA strand
            - transcript_id: ID of the transcript
            - is_cds: Boolean indicating if feature is CDS
        """
        try:
            fields = line.strip().split('\t')
            if len(fields) < 9:
                return None

            attributes = dict(attr.strip().split(' ', 1) for attr in
                          fields[8].rstrip(';').split('; ') if ' ' in attr)
            transcript_id = attributes.get('transcript_id', '').strip('"')
            feature_type = fields[2]

            return {
                'chromosome': fields[0],
                'feature': feature_type,
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'transcript_id': transcript_id.split('.')[0],
                'is_cds': feature_type == 'CDS'
            }
        except Exception as e:
            logging.error(f"Error parsing GTF line: {e}")
            return None