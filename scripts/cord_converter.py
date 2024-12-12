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
									'exons': []
								}
							self.transcripts[transcript_id]['exons'].append((int(fields[3]), int(fields[4])))

		# Sort exons by coordinates
		for transcript in self.transcripts.values():
			transcript['exons'].sort()

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
		Converts genomic coordinates to transcript coordinates.
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

		transcript_pos = 0
		for start, end in transcript['exons']:
			if start <= pos <= end:
				if transcript['strand'] == '+':
					transcript_pos += pos - start + 1
				else:
					transcript_pos += end - pos + 1
				return transcript_pos
			transcript_pos += end - start + 1

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
			chrom, pos = row[0], row[1]
			transcript_field = row[11]
			transcript_id = transcript_field.split('|')[0]
			if transcript_id not in self.converter.transcripts:
				print(f"Warning: Transcript {transcript_id} not found in the GTF.")
				continue
			transcript_coord = self.converter.genomic_to_transcript(chrom, int(pos), transcript_id)
			print(f"Genomic: {chrom}:{pos}, Transcript: {transcript_id}:{transcript_coord}")


# Example usage
pipeline = Pipeline("data/sorted.v4.bed", "data/HGMD_resort.vcf", "data/combined_uorf_sorted.v4.gtf.gz")
pipeline.run()
