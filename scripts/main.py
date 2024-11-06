from data_loader import GTFFile, VCFFile, BEDFile
from data_processor import Bedtools, DataProcessor


class Pipeline:
	def __init__(self, gtf_file_path, vcf_file_path, bed_file_path):
		self.gtf = GTFFile(gtf_file_path)
		self.vcf = VCFFile(vcf_file_path)
		self.bed = BEDFile(bed_file_path)
		self.intersected_bed = ''
		self.processor = DataProcessor()

	def run(self):
		self._initial_process_gtf_file()
		self._process_vcf_file()
		self._process_bed_file()
		self._intersect_vcf_and_bed()
		self._process_data()
		# print(self.gtf.uorf_df)

	def _initial_process_gtf_file(self):
		self.gtf.process_gtf_file()

	def _process_vcf_file(self):
		return self.vcf.process_vcf_file()

	def _process_bed_file(self):
		return self.bed.process_bed_file()

	def _intersect_vcf_and_bed(self):
		self.intersected_bed = Bedtools.intersect(self.vcf.file_path, self.bed.file_path)

	def _process_data(self):
		intersection_file_path = self.intersected_bed
		bed_4col_info_cols = self.bed.bed_4col_info_cols
		gene_transcript_records = self.gtf.gene_transcript
		source_gtf = self.gtf.source_gtf
		self.processor.process_data(intersection_file_path, bed_4col_info_cols, gene_transcript_records, source_gtf)



if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Extract gene_id and transcript_id from a GTF file, read the header and table of a VCF file, and intersect with a BED file.')
	parser.add_argument('gtf_file', type=str, help='Path to the GTF file.')
	parser.add_argument('vcf_file', type=str, help='Path to the VCF file.')
	parser.add_argument('bed_file', type=str, help='Path to the BED file.')
	args = parser.parse_args()

	pipeline = Pipeline(args.gtf_file, args.vcf_file, args.bed_file)
	result = pipeline.run()
