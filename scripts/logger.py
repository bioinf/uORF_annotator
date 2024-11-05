class Logger:
	@staticmethod
	def log_gene_transcript_records(num_records):
		print(f"Number of gene-transcript records from GTF: {num_records}")

	@staticmethod
	def log_processed_bed_records(num_records):
		print(f"Number of processed intervals in BED file: {num_records}")

	@staticmethod
	def log_num_variants_in_intersection(num_variants):
		print(f"Number of variants in intersection: {num_variants}")

	@staticmethod
	def log_num_records_in_table_after_annotation_processing_1(shape):
		print(f"Number of rows and columns in table after annotation processing #1: {shape}")

	@staticmethod
	def log_num_exons_after_gtf_processing(n_rows):
		print(f"Number of exons in BED file after GTF processing: {n_rows}")
	
	@staticmethod
	def log_num_records_in_table_after_annotation_processing_2(shape):
		print(f"Number of rows and columns in table after uORF extraction: {shape}")

	@staticmethod
	def log_num_records_in_table_after_annotation_processing_3(shape):
		print(f"Number of rows and columns in table after interorf_single processing: {shape}")
