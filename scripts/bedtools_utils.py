# bedtools_utils.py

import logging
from typing import Optional, Dict
from pybedtools import BedTool
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class BedToolsWrapper:
    """Wrapper for BedTools operations."""
    
    @staticmethod
    def intersect(vcf_file: str, bed_file: str) -> pd.DataFrame:
        """
        Intersect VCF and BED files using pybedtools.
        
        Args:
            vcf_file: Path to VCF file
            bed_file: Path to BED file
            
        Returns:
            DataFrame with intersected regions
        """
        logger.info("Creating BedTool objects...")
        vcf = BedTool(vcf_file)
        bed = BedTool(bed_file)

        logger.info("Performing intersection...")
        intersection = vcf.intersect(bed, wa=True, wb=True)
        
        logger.info("Converting intersection to DataFrame...")
        columns = [
            'chrom', 'pos', 'id', 'ref', 'alt', 'name', 'score', 'info',
            'uorf_chrom', 'uorf_start', 'uorf_end', 'uorf_info', 'uorf_a',
            'uorf_strand', 'bed_start', 'bed_end', 'uorf_b', 'uorf_c', 'exon_sizes', 'exon_starts'
        ]
        return pd.read_csv(intersection.fn, sep='\t', index_col=False, header=None, names=columns)

    @staticmethod
    def get_fasta(bed_file: str, fasta_file: str, strand: bool = True, name: bool = True) -> Dict[str, str]:
        """
        Extract sequences from FASTA file based on BED coordinates.
        
        Args:
            bed_file: Path to BED file
            fasta_file: Path to FASTA file
            strand: Consider strand information
            name: Use names in BED file
            
        Returns:
            Dictionary mapping transcript IDs to their sequences
        """
        logger.info(f"Extracting sequences from {fasta_file} using coordinates from {bed_file}")
        
        try:
            # Create BedTool object and extract sequences
            bed = BedTool(bed_file)
            sequences = bed.sequence(fi=fasta_file, s=strand, name=name)
            
            # Parse FASTA output into dictionary
            result = {}
            current_id = None
            current_seq = []
            
            with open(sequences.seqfn) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            result[current_id] = ''.join(current_seq)
                            current_seq = []
                        # Parse transcript ID from BED name field
                        current_id = line[1:].split('|')[0].split('.')[0]  # Get NM_XXXXX part
                    else:
                        current_seq.append(line)
                        
            # Add last sequence
            if current_id:
                result[current_id] = ''.join(current_seq)
                
            logger.info(f"Successfully extracted {len(result)} sequences")
            
            # Добавляем дополнительное логирование
            for transcript_id, seq in result.items():
                logger.debug(f"Extracted sequence for {transcript_id}: length={len(seq)}, first 50bp={seq[:50]}")
                
            return result
            
        except Exception as e:
            logger.error(f"Error extracting sequences: {str(e)}")
            raise

    @staticmethod
    def check_file_format(file_path: str, file_type: str) -> bool:
        """
        Validate file format.
        
        Args:
            file_path: Path to file
            file_type: Expected file type ('bed', 'vcf', or 'fasta')
            
        Returns:
            True if file format is valid, False otherwise
        """
        try:
            if file_type == 'bed':
                BedTool(file_path)
            elif file_type == 'vcf':
                with open(file_path) as f:
                    for line in f:
                        if not line.startswith('#'):
                            fields = line.strip().split('\t')
                            if len(fields) < 8:
                                return False
                            break
            elif file_type == 'fasta':
                with open(file_path) as f:
                    first_line = f.readline().strip()
                    if not first_line.startswith('>'):
                        return False
            return True
            
        except Exception as e:
            logger.error(f"Error validating {file_type} file: {str(e)}")
            return False