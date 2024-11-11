from data_processor import DataProcessor

import argparse


if __name__ == "__main__":

	"python scripts/main.py data/combined_uorf.v4.gtf.gz data/HGMD_resort.vcf data/sorted.v4.bed"

	parser = argparse.ArgumentParser(description='Extract gene_id and transcript_id from a GTF file, read the header and table of a VCF file, and intersect with a BED file.')
	parser.add_argument('gtf_file', type=str, help='Path to the GTF file.')
	parser.add_argument('vcf_file', type=str, help='Path to the VCF file.')
	parser.add_argument('bed_file', type=str, help='Path to the BED file.')
	parser.add_argument('fasta_file', type=str, help='Path to the FASTA file.')
	args = parser.parse_args()

	pipeline = DataProcessor(args.gtf_file, args.vcf_file, args.bed_file, args.fasta_file)
	result = pipeline.process_data()
