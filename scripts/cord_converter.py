import pysam
from pybedtools import BedTool
import pandas as pd

class CoordinateConverter:
    def __init__(self, gtf_file):
        """
        Initializes the converter using a GTF file.
        :param gtf_file: path to the GTF file with gene annotations.
        """
        self.gtf_file = gtf_file
        self.transcripts = {}
        self._parse_gtf()

    def _parse_gtf(self):
        """
        Parses the GTF file to extract information about transcripts and exons.
        """
        with pysam.TabixFile(self.gtf_file) as gtf:
            for chrom in gtf.contigs:
                for line in gtf.fetch(chrom):
                    fields = line.split('\t')
                    if fields[2] == 'exon':
                        attributes = self._parse_attributes(fields[8])
                        transcript_id = attributes.get('transcript_id')
                        if transcript_id:
                            if transcript_id not in self.transcripts:
                                self.transcripts[transcript_id] = {
                                    'chrom': fields[0],
                                    'strand': fields[6],
                                    'exons': [],
                                    'cumulative_lengths': []
                                }
                            self.transcripts[transcript_id]['exons'].append((int(fields[3]), int(fields[4])))

        # Sort exons and calculate cumulative lengths
        for transcript_id, transcript in self.transcripts.items():
            # Sort exons based on genomic coordinates
            transcript['exons'].sort(key=lambda x: x[0])
            
            # Calculate total transcript length
            total_length = sum(end - start + 1 for start, end in transcript['exons'])
            transcript['total_length'] = total_length
            
            # Calculate cumulative lengths
            cumulative_length = 0
            cumulative_lengths = []
            
            if transcript['strand'] == '+':
                for start, end in transcript['exons']:
                    exon_length = end - start + 1
                    cumulative_lengths.append((cumulative_length, cumulative_length + exon_length))
                    cumulative_length += exon_length
            else:
                # For negative strand, we need to calculate cumulative lengths in reverse
                for start, end in reversed(transcript['exons']):
                    exon_length = end - start + 1
                    cumulative_length_start = total_length - cumulative_length - exon_length
                    cumulative_lengths.append((cumulative_length_start, cumulative_length_start + exon_length))
                    cumulative_length += exon_length
                
            transcript['cumulative_lengths'] = cumulative_lengths

    def _parse_attributes(self, attribute_string):
        """
        Parses the GTF attribute string and returns a dictionary.
        :param attribute_string: GTF attribute string
        :return: dictionary with attributes
        """
        attributes = {}
        for attr in attribute_string.strip().split(';'):
            if attr.strip():
                parts = attr.strip().split(' ', 1)
                if len(parts) == 2:
                    key, value = parts
                    attributes[key] = value.strip('"')
        return attributes

    def genomic_to_transcript(self, chrom, pos, transcript_id):
        """
        Converts genomic coordinates to transcript coordinates, considering exon-intron structure.
        :param chrom: chromosome
        :param pos: genomic position (1-based)
        :param transcript_id: transcript identifier
        :return: transcript coordinate or None if the position is not within exons
        """
        if transcript_id not in self.transcripts:
            raise ValueError(f"Transcript {transcript_id} not found in the annotation.")
            
        transcript = self.transcripts[transcript_id]
        if transcript['chrom'] != chrom:
            return None

        # Convert position to 0-based for calculations
        pos_0based = pos - 1

        if transcript['strand'] == '+':
            # For positive strand
            for (exon_start, exon_end), (cum_start, cum_end) in zip(
                transcript['exons'],
                transcript['cumulative_lengths']
            ):
                if exon_start - 1 <= pos_0based <= exon_end - 1:
                    relative_pos = pos_0based - (exon_start - 1)
                    return cum_start + relative_pos + 1  # Convert back to 1-based
        else:
            # For negative strand
            for (exon_start, exon_end), (cum_start, cum_end) in zip(
                reversed(transcript['exons']),
                transcript['cumulative_lengths']
            ):
                if exon_start - 1 <= pos_0based <= exon_end - 1:
                    relative_pos = exon_end - 1 - pos_0based
                    return cum_start + relative_pos + 1  # Convert back to 1-based
                    
        return None

    def transcript_to_genomic(self, transcript_id, transcript_pos):
        """
        Converts transcript coordinates to genomic coordinates.
        :param transcript_id: transcript identifier
        :param transcript_pos: position in transcript coordinates (1-based)
        :return: tuple (chrom, genomic_pos) or None if position is invalid
        """
        if transcript_id not in self.transcripts:
            raise ValueError(f"Transcript {transcript_id} not found in the annotation.")
            
        transcript = self.transcripts[transcript_id]
        
        # Check if position is within transcript bounds
        if transcript_pos < 1 or transcript_pos > transcript['total_length']:
            return None
        
        # Convert to 0-based for calculations
        transcript_pos_0based = transcript_pos - 1
        
        if transcript['strand'] == '+':
            # For positive strand, scan exons in forward order
            for (exon_start, exon_end), (cum_start, cum_end) in zip(
                transcript['exons'],
                transcript['cumulative_lengths']
            ):
                if cum_start <= transcript_pos_0based < cum_end:
                    relative_pos = transcript_pos_0based - cum_start
                    genomic_pos = exon_start + relative_pos
                    return (transcript['chrom'], genomic_pos)
        else:
            # For negative strand, scan exons in reverse order
            for (exon_start, exon_end), (cum_start, cum_end) in zip(
                reversed(transcript['exons']),
                transcript['cumulative_lengths']
            ):
                if cum_start <= transcript_pos_0based < cum_end:
                    relative_pos = cum_end - transcript_pos_0based - 1
                    genomic_pos = exon_end - relative_pos
                    return (transcript['chrom'], genomic_pos)
                
        return None

class Bedtools:
    @staticmethod
    def intersect(file1: str, file2: str) -> pd.DataFrame:
        """
        Performs an intersection between two BED files.
        :param file1: path to the first BED file
        :param file2: path to the second BED file
        :return: a DataFrame containing the intersection
        """
        bed_file1 = BedTool(file1)
        bed_file2 = BedTool(file2)
        intersected_bed = bed_file1.intersect(bed_file2, wo=True)
        intersection_df = intersected_bed.to_dataframe(header=None)
        intersection_df.columns = range(len(intersection_df.columns))
        return intersection_df

    @staticmethod
    def getfasta(fasta_file: str, bed_file: str, strand: bool = True) -> str:
        """
        Extracts sequences from a FASTA file based on a BED file.
        :param fasta_file: path to the FASTA file
        :param bed_file: path to the BED file
        :param strand: whether to consider the strand information
        :return: path to the output FASTA file
        """
        fasta = BedTool(fasta_file)
        bed = BedTool(bed_file)
        getfasta = bed.sequence(fi=fasta, s=strand)
        return getfasta.seqfn

class Pipeline:
    def __init__(self, bed_file: str, vcf_file: str, gtf_file: str):
        """
        Initializes the pipeline with the necessary files.
        :param bed_file: path to the BED file
        :param vcf_file: path to the VCF file
        :param gtf_file: path to the GTF file
        """
        self.bed_file = bed_file
        self.vcf_file = vcf_file
        self.gtf_file = gtf_file
        self.converter = CoordinateConverter(gtf_file)

    def run(self):
        """
        Runs the pipeline logic step by step.
        """
        intersection = Bedtools.intersect(self.vcf_file, self.bed_file)

        for _, row in intersection.iterrows():
            chrom, pos = row[0], int(row[1])
            transcript_field = row[11]
            transcript_id = transcript_field.split('|')[0]
            
            if transcript_id not in self.converter.transcripts:
                print(f"Warning: Transcript {transcript_id} not found in the GTF.")
                continue
                
            # Convert genomic to transcript coordinates
            transcript_coord = self.converter.genomic_to_transcript(chrom, pos, transcript_id)
            if transcript_coord is None:
                print(f"Warning: Position {chrom}:{pos} is not within any exon of transcript {transcript_id}")
                continue
            
            print(f"Genomic: {chrom}:{pos}, Transcript: {transcript_id}:{transcript_coord}")
            
            # Verify conversion by converting back to genomic coordinates
            back_conversion = self.converter.transcript_to_genomic(transcript_id, transcript_coord)
            if back_conversion:
                chrom_back, pos_back = back_conversion
                print(f"Verification - converting back: Transcript: {transcript_id}:{transcript_coord} -> Genomic: {chrom_back}:{pos_back}")
                if pos != pos_back:
                    print(f"Warning: Coordinate conversion mismatch! Original: {pos}, After conversion: {pos_back}")

# Example usage
pipeline = Pipeline("data/sorted.v4.bed", "data/HGMD_resort.vcf", "data/combined_uorf_sorted.v4.gtf.gz")
pipeline.run()