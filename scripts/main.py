from file_processor import GTFFile, VCFFile, BEDFile
from annotator import BedtoolsIntersect, AnnotationProcessor


class Pipeline:
	def __init__(self, gtf_file_path, vcf_file_path, bed_file_path):
		self.gtf_file = GTFFile(gtf_file_path)
		self.vcf_file = VCFFile(vcf_file_path)
		self.bed_file = BEDFile(bed_file_path)
		self.intersector = BedtoolsIntersect("/tmp")
		self.processor = AnnotationProcessor()

	def run(self):
		gene_transcript_records = self._get_gene_transcript_records()
		vcf_header = self._read_vcf_header()
		source_bed_lines = self._process_bed_file()
		intersection_file = self._intersect_vcf_and_bed()
		bed_4col_info_cols = self._get_bed_4col_info_cols()
		data_annotated = self._process_annotation(intersection_file.name, bed_4col_info_cols, gene_transcript_records)
		exons_file = self._process_gtf_file()

	def _get_gene_transcript_records(self):
		return self.gtf_file.get_gene_transcript_records()

	def _read_vcf_header(self):
		return self.vcf_file.read_header()

	def _process_bed_file(self):
		return self.bed_file.process_bed_file()

	def _intersect_vcf_and_bed(self):
		return self.intersector.intersect(self.vcf_file.file_path, self.bed_file.file_path)

	def _get_bed_4col_info_cols(self):
		return self.bed_file.bed_4col_info_cols

	def _process_annotation(self, intersection_file_path, bed_4col_info_cols, gene_transcript_records):
		return self.processor.process_data(intersection_file_path, bed_4col_info_cols, gene_transcript_records)

	def _process_gtf_file(self):
		return self.gtf_file.process_gtf_file()


if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Extract gene_id and transcript_id from a GTF file, read the header and table of a VCF file, and intersect with a BED file.')
	parser.add_argument('gtf_file', type=str, help='Path to the GTF file.')
	parser.add_argument('vcf_file', type=str, help='Path to the VCF file.')
	parser.add_argument('bed_file', type=str, help='Path to the BED file.')
	args = parser.parse_args()

	pipeline = Pipeline(args.gtf_file, args.vcf_file, args.bed_file)
	result = pipeline.run()
