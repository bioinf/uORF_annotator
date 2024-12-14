import argparse
import pandas as pd
import logging
from dataclasses import dataclass
from typing import List, Tuple, Dict
from pybedtools import BedTool
import gzip
import csv


# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class Exon:
    start: int
    length: int
    genome_start: int
    genome_end: int


class Transcript:
    def __init__(self, transcript_id: str, chromosome: str, strand: str, exons: List[Exon]):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = sorted(exons, key=lambda x: x.genome_start)
        self._build_coordinate_maps()

    def _build_coordinate_maps(self):
        """Build mappings between genomic and transcript coordinates."""
        # logger.debug(f"Building coordinate maps for transcript {self.transcript_id}")
        self.genome_to_transcript = {}
        self.transcript_to_genome = {}

        current_transcript_pos = 0
        for exon in self.exons:
            for genome_pos in range(exon.genome_start, exon.genome_end + 1):
                if self.strand == '+':
                    self.genome_to_transcript[genome_pos] = current_transcript_pos
                    self.transcript_to_genome[current_transcript_pos] = genome_pos
                else:
                    self.genome_to_transcript[genome_pos] = current_transcript_pos
                    self.transcript_to_genome[current_transcript_pos] = genome_pos
                current_transcript_pos += 1


class BedToolsWrapper:
    @staticmethod
    def intersect(vcf_file: str, bed_file: str) -> pd.DataFrame:
        """Intersect VCF and BED files using pybedtools."""
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
        df = pd.read_csv(intersection.fn, sep='\t', index_col=False, header=None, names=columns)
        return df


class CoordinateConverter:
    def __init__(self, bed_file: str, gtf_file: str):
        self.bed_file = bed_file
        self.gtf_file = gtf_file
        self.transcripts: Dict[str, Transcript] = {}
        self._parse_files()

    def _parse_gtf_line(self, line: str) -> dict:
        """Parse a single GTF line and return relevant information."""
        try:
            fields = line.strip().split('\t')
            if len(fields) < 9:
                return None

            attributes = dict(attr.strip().split(' ', 1) for attr in
                              fields[8].rstrip(';').split('; ') if ' ' in attr)
            transcript_id = attributes.get('transcript_id', '').strip('"')

            return {
                'chromosome': fields[0],
                'feature': fields[2],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'transcript_id': transcript_id.split('.')[0]
            }
        except Exception as e:
            logger.error(f"Error parsing GTF line: {e}")
            return None

    def _parse_files(self):
        """Parse BED and GTF files to build transcript information."""
        logger.info("Starting to parse input files...")

        logger.info("Reading BED file...")
        bed_df = pd.read_csv(self.bed_file, sep='\t', header=None)

        exons_by_transcript = {}

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
                if not parsed or parsed['feature'] != 'exon':
                    continue

                transcript_id = parsed['transcript_id']
                if not transcript_id:
                    continue

                if transcript_id not in exons_by_transcript:
                    exons_by_transcript[transcript_id] = []

                exon = Exon(
                    start=len(exons_by_transcript[transcript_id]),
                    length=parsed['end'] - parsed['start'] + 1,
                    genome_start=parsed['start'],
                    genome_end=parsed['end']
                )
                exons_by_transcript[transcript_id].append(exon)

        logger.info("Creating Transcript objects...")
        for _, bed_row in bed_df.iterrows():
            transcript_id = bed_row[3].split('|')[0].split('.')[0]
            if transcript_id in exons_by_transcript:
                self.transcripts[transcript_id] = Transcript(
                    transcript_id=transcript_id,
                    chromosome=bed_row[0],
                    strand=bed_row[5],
                    exons=exons_by_transcript[transcript_id]
                )

        logger.info(f"Created {len(self.transcripts)} transcript objects")

    def genome_to_transcript_pos(self, transcript_id: str, genome_pos: int) -> int:
        """Convert genomic position to transcript position."""
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            logger.warning(f"Transcript {transcript_id} not found")
            return -1

        return transcript.genome_to_transcript.get(genome_pos, -1)

    def transcript_to_genome_pos(self, transcript_id: str, transcript_pos: int) -> int:
        """Convert transcript position to genomic position."""
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            logger.warning(f"Transcript {transcript_id} not found")
            return -1

        return transcript.transcript_to_genome.get(transcript_pos, -1)


class Pipeline:
    def __init__(self, bed_file: str, vcf_file: str, gtf_file: str):
        self.bed_file = bed_file
        self.vcf_file = vcf_file
        self.gtf_file = gtf_file
        self.bedtools = BedToolsWrapper()
        logger.info("Initializing CoordinateConverter...")
        self.converter = CoordinateConverter(bed_file, gtf_file)

    def check_position_in_exons(self, transcript_id: str, genome_pos: int) -> bool:
        """Check if genomic position falls within any exon of the transcript."""
        transcript = self.converter.transcripts.get(transcript_id)
        if not transcript:
            return False
        
        for exon in transcript.exons:
            if exon.genome_start <= genome_pos <= exon.genome_end:
                return True
        return False

    def get_nearest_exon_boundaries(self, transcript_id: str, genome_pos: int) -> Tuple[int, int]:
        """Find the nearest exon boundaries for a given genomic position."""
        transcript = self.converter.transcripts.get(transcript_id)
        if not transcript:
            return None, None
        
        prev_end = None
        next_start = None
        
        for exon in transcript.exons:
            if genome_pos < exon.genome_start:
                next_start = exon.genome_start
                break
            elif exon.genome_start <= genome_pos <= exon.genome_end:
                return exon.genome_start, exon.genome_end
            prev_end = exon.genome_end
            
        return prev_end, next_start

    def process_variants(self) -> pd.DataFrame:
        """Process VCF variants and convert coordinates."""
        logger.info("Intersecting VCF with BED...")
        intersected = self.bedtools.intersect(self.vcf_file, self.bed_file)

        logger.info("Converting coordinates...")
        results = []
        total_rows = len(intersected)
        failed_conversions = 0

        for i in range(total_rows):
            if (i >= 0 and i % 1000 == 0) or (i == total_rows - 1):
                logger.info(f"Processed {i}/{total_rows} variants")

            row = intersected.iloc[i]
            vcf_pos = row['pos']
            transcript_id = row['uorf_info'].split('|')[0].split('.')[0]

            # Check that transcript_id exists in self.transcripts dict
            if transcript_id not in self.converter.transcripts:
                logger.warning(f"Transcript {transcript_id} not found")
                continue

            transcript_pos = self.converter.genome_to_transcript_pos(transcript_id, vcf_pos)
            
            if transcript_pos == -1:
                failed_conversions += 1
                in_exon = self.check_position_in_exons(transcript_id, vcf_pos)
                prev_end, next_start = self.get_nearest_exon_boundaries(transcript_id, vcf_pos)
                
                logger.warning(
                    f"Failed to convert position {vcf_pos} for transcript {transcript_id}:\n"
                    f"  Position in exon: {in_exon}\n"
                    f"  Nearest exon boundaries: {prev_end} - {next_start}\n"
                    f"  VCF record: {row['chrom']}:{vcf_pos} {row['ref']} {row['alt']}\n"
                    f"  Strand: {row['uorf_strand']}"
                )
            
            # Convert transcript position back to genome coordinates for validation
            recovered_genome_pos = self.converter.transcript_to_genome_pos(transcript_id, transcript_pos)

            results.append({
                'Chromosome': row['chrom'],
                'Original_Genome_Position': vcf_pos,
                'Ref_Allele': row['ref'],
                'Alt_Allele': row['alt'],
                'Strand': row['uorf_strand'],
                'Transcript_ID': transcript_id,
                'Transcript_Position': transcript_pos,
                'Recovered_Genome_Position': recovered_genome_pos,
                'Coordinates_Match': vcf_pos == recovered_genome_pos,
                'In_Exon': self.check_position_in_exons(transcript_id, vcf_pos),
                'Nearest_Exon_Start': self.get_nearest_exon_boundaries(transcript_id, vcf_pos)[1],
                'Nearest_Exon_End': self.get_nearest_exon_boundaries(transcript_id, vcf_pos)[0]
            })

        if failed_conversions > 0:
            logger.warning(f"Total failed conversions: {failed_conversions} out of {total_rows} "
                         f"({(failed_conversions/total_rows)*100:.2f}%)")

        return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description='Convert genomic coordinates to transcript coordinates and vice versa')
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--output', required=True, help='Path to output file')
    # parser.add_argument('--debug', action='store_true', help='Enable debug logging')

    args = parser.parse_args()

    # if args.debug:
    #     logger.setLevel(logging.DEBUG)

    try:
        logger.info("Starting pipeline...")
        pipeline = Pipeline(args.bed, args.vcf, args.gtf)
        results = pipeline.process_variants()
        results.to_csv(args.output, sep='\t', index=False)
        logger.info(f"Results saved to {args.output}")
        
        # Log coordinate matching statistics
        total_variants = len(results)
        matching_coords = results['Coordinates_Match'].sum()
        logger.info(f"Coordinate conversion validation:")
        logger.info(f"Total variants processed: {total_variants}")
        logger.info(f"Variants with matching coordinates: {matching_coords} ({(matching_coords/total_variants)*100:.2f}%)")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
        raise
    finally:
        logger.info("Pipeline finished successfully.")


if __name__ == '__main__':
    main()