from logger import Logger
from temporary_file_manager import TemporaryFileManager


import pandas as pd
import subprocess
import re


class BedtoolsIntersect:
	def __init__(self, tmp_dir):
		self.tmp_dir = tmp_dir

	def intersect(self, file1, file2):
		tmp_out = TemporaryFileManager.create('.bed', self.tmp_dir)
		command = f"bedtools intersect -wo -a {file1} -b {file2}"
		with open(tmp_out.name, 'w') as w:
			subprocess.run(command, shell=True, stdout=w)
		num_lines = int(subprocess.check_output(f"wc -l {tmp_out.name}", shell=True).split()[0])
		Logger.log_num_variants_in_intersection(num_lines)

		tmp_out.register_at_exit()

		return tmp_out

class AnnotationProcessor:
	def __init__(self):
		self.annotation_data = None

	def process_data(self, intersection_file_path, bed_4col_info_cols, gene_transcript_records):
		self._load_data(intersection_file_path)
		self._extract_exons()
		self._extract_bed_anno(bed_4col_info_cols)
		self._extract_positions()
		self._extract_names(gene_transcript_records)
		Logger.log_num_records_in_table_after_annotation_processing_1(self.annotation_data.shape)
		return(self.annotation_data)

	def _load_data(self, intersection_file_path):
		self.annotation_data = pd.read_table(intersection_file_path, header=None)
		self.annotation_data = self.annotation_data.dropna()
		self.annotation_data = self.annotation_data.loc[:, [0, 1, 3, 4, 7, 9, 10, 11, 13, 17, 18, 19]]
		self.annotation_data = self.annotation_data.rename({0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT',
													7: 'INFO', 9: 'orf_start', 10: 'orf_end',
													11: 'bed_anno', 13:'strand', 17:'n_exons', 
													18: 'exons_sizes', 19: 'exons_starts'}, axis=1)

	def _extract_exons(self):
		self.annotation_data['exons_sizes'] = self.annotation_data['exons_sizes'].apply(lambda x: list(map(int, x.rstrip(',').split(','))))
		self.annotation_data['exons_starts'] = self.annotation_data['exons_starts'].apply(lambda x: list(map(int, x.rstrip(',').split(','))))

	def _extract_bed_anno(self, bed_4col_info_cols):
		self.annotation_data = pd.concat([self.annotation_data, self.annotation_data['bed_anno'].str.split('|', expand=True)], axis=1)
		self.annotation_data.columns = list(self.annotation_data.columns[:-4]) + bed_4col_info_cols

	def _extract_positions(self):
		self.annotation_data[['POS', 'orf_start', 'orf_end']] = self.annotation_data[['POS', 'orf_start', 'orf_end']].astype(int)
		self.annotation_data['dist_from_orf_to_snp'] = self.annotation_data.apply(lambda x: x['POS'] - x['orf_start'] - 1 if x['strand'] == '+' else x['orf_end'] - x['POS'], axis=1).astype(int)

	def _extract_names(self, gene_transcript_records):
		self.annotation_data['name'] = self.annotation_data.apply(lambda x: f"{x['#CHROM']}:{x['orf_start']}-{x['orf_end']}({x['strand']})", axis=1)
		self.annotation_data['transcript'] = self.annotation_data.apply(lambda x: re.findall('[NMENSTX]+_*\d+', x['bed_anno'])[0], axis=1)
		self.annotation_data['gene_name'] = self.annotation_data.apply(lambda x: gene_transcript_records.get(x['transcript'], ''), axis=1)
		self.annotation_data['name_and_trans'] = self.annotation_data.apply(lambda x: f"{x['name']}/{x['transcript']}", axis=1)
