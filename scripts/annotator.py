from logger import Logger
from temporary_file_manager import TemporaryFileManager


import pandas as pd
import subprocess
import re

import hashlib


class Bedtools:
	@staticmethod
	def intersect(file1, file2, tmp_dir='/tmp'):
		bed_file = TemporaryFileManager.create('.bed', tmp_dir)
		command = f"bedtools intersect -wo -a {file1} -b {file2}"
		with open(bed_file.name, 'w') as w:
			subprocess.run(command, shell=True, stdout=w)
		num_lines = int(subprocess.check_output(f"wc -l {bed_file.name}", shell=True).split()[0])

		Logger.log_num_variants_in_intersection(num_lines)

		bed_file.register_at_exit()

		return bed_file.name


class DataProcessor:
	def __init__(self):
		self.data = None
		self.first_cds_df = None
		self.uorf_data = None
		self.interorf_single = None

	def process_data(self, intersection_file_path, bed_4col_info_cols, gene_transcript_records, source_gtf):
		self._load_data(intersection_file_path)
		self._extract_exon_starts_and_sizes()
		self._extract_predefined_bed_anno(bed_4col_info_cols)
		self._extract_distancies_from_orf_to_snp()
		self._extract_names(gene_transcript_records)
		self._logger_1()
		self._get_first_cds(source_gtf)
		self._extract_uorf_data()
		self._logger_2()
		self._extract_interorf_data()
		self._logger_3()
		md5sum_uorf = hashlib.md5(self.uorf_data.to_csv(index=False).encode()).hexdigest()
		md5sum_interorf = hashlib.md5(self.interorf_single.to_csv(index=False).encode()).hexdigest()
		print("md5sum:", md5sum_uorf)
		print("md5sum:", md5sum_interorf)

	def _logger_1(self):
		Logger.log_num_records_in_table_after_annotation_processing_1(self.data.shape)

	def _logger_2(self):
		Logger.log_num_records_in_table_after_annotation_processing_2(self.uorf_data.shape)

	def _logger_3(self):
		Logger.log_num_records_in_table_after_annotation_processing_3(self.interorf_single.shape)

	def _load_data(self, intersection_file_path):
		self.data = pd.read_table(intersection_file_path, header=None)
		self.data = self.data.dropna()
		self.data = self.data.loc[:, [0, 1, 3, 4, 7, 9, 10, 11, 13, 17, 18, 19]]
		self.data = self.data.rename({0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT',
													7: 'INFO', 9: 'orf_start', 10: 'orf_end',
													11: 'bed_anno', 13:'strand', 17:'n_exons', 
													18: 'exons_sizes', 19: 'exons_starts'}, axis=1)

	def _extract_exon_starts_and_sizes(self):
		self.data['exons_starts'] = self.data['exons_starts'].apply(lambda x: list(map(int, x.rstrip(',').split(','))))
		self.data['exons_sizes'] = self.data['exons_sizes'].apply(lambda x: list(map(int, x.rstrip(',').split(','))))

	def _extract_predefined_bed_anno(self, bed_4col_info_cols):
		self.data = pd.concat([self.data, self.data['bed_anno'].str.split('|', expand=True)], axis=1)
		self.data.columns = list(self.data.columns[:-4]) + bed_4col_info_cols

	def _extract_distancies_from_orf_to_snp(self):
		self.data[['POS', 'orf_start', 'orf_end']] = self.data[['POS', 'orf_start', 'orf_end']].astype(int)
		self.data['dist_from_orf_to_snp'] = self.data.apply(lambda x: x['POS'] - x['orf_start'] - 1 if x['strand'] == '+' else x['orf_end'] - x['POS'], axis=1).astype(int)

	def _extract_names(self, gene_transcript_records):
		self.data['name'] = self.data.apply(lambda x: f"{x['#CHROM']}:{x['orf_start']}-{x['orf_end']}({x['strand']})", axis=1)
		self.data['transcript'] = self.data.apply(lambda x: re.findall('[NMENSTX]+_*\d+', x['bed_anno'])[0], axis=1)
		self.data['gene_name'] = self.data.apply(lambda x: gene_transcript_records.get(x['transcript'], ''), axis=1)
		self.data['name_and_trans'] = self.data.apply(lambda x: f"{x['name']}/{x['transcript']}", axis=1)
	
	def _get_first_cds(self, source_gtf):
		self.first_cds_df = source_gtf.loc[source_gtf['type'] == 'CDS'].sort_values(['transcript', 'strand', 'start']).groupby('transcript').apply(lambda x: x.iloc[0] if x['strand'].iloc[0] == '+' else x.iloc[-1]).reset_index(drop=True)
		self.first_cds_df.index = self.first_cds_df['transcript']

	def _get_true_cds_start(self, tr_id):
		return self.first_cds_df.loc[tr_id, 'start'] if tr_id in self.first_cds_df.index and self.first_cds_df.loc[tr_id, 'strand'] == '+' else self.first_cds_df.loc[tr_id, 'end'] if tr_id in self.first_cds_df.index else None

	def _extract_uorf_data(self):
		self.uorf_data = self.data.loc[:, ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_starts', 'name', 'name_and_trans']].copy()
		self.uorf_data.columns = ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_norm_starts', 'id', 'name_and_trans']
		self.uorf_data = self.uorf_data.applymap(lambda x: str(x) if isinstance(x, list) else x).drop_duplicates()
		self.uorf_data['exons_norm_starts'] = self.uorf_data['exons_norm_starts'].apply(eval)
		self.uorf_data['exons_sizes'] = self.uorf_data['exons_sizes'].apply(eval)
		self.uorf_data[['exons_starts', 'exons_norm_starts']] = self.uorf_data.apply(self._process_exons, axis=1, result_type='expand')
		self.uorf_data['exons_sizes'] = self.uorf_data.apply(lambda row: row['exons_sizes'] if row['strand'] == '+' else row['exons_sizes'][::-1], axis=1)
		self.uorf_data.index = self.uorf_data['id']
		self.uorf_data['cds_start'] = self.uorf_data['transcript'].apply(self._get_true_cds_start)

	def _process_exons(self, row):
		orf_len = row['orf_end'] - row['orf_start']
		exons_norm_starts = row['exons_norm_starts']
		exons_sizes = row['exons_sizes']
		if row['strand'] == '+':
			return [row['orf_start'] + eStart for eStart in exons_norm_starts], exons_norm_starts
		else:
			return [row['orf_start'] + eStart + eSize for eStart, eSize in zip(exons_norm_starts, exons_sizes)][::-1], [orf_len - eStart - eSize for eStart, eSize in zip(exons_norm_starts, exons_sizes)][::-1]

	def _extract_interorf_data(self):
		mask = ((self.uorf_data['strand'] == '+') & (self.uorf_data['orf_end'] < self.uorf_data['cds_start'])) | (self.uorf_data['cds_start'] < self.uorf_data['orf_start'])
		interorf_data = self.uorf_data[mask].copy()
		interorf_data.loc[interorf_data['strand'] == '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] == '+', ['orf_end', 'cds_start']].values
		interorf_data.loc[interorf_data['strand'] != '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] != '+', ['cds_start', 'orf_start']].values
		self.interorf_single = interorf_data[['#CHROM', 'start', 'end', 'name_and_trans']].rename(columns={'#CHROM': 'chr'}).dropna()
		self.interorf_single[['start', 'end']] = self.interorf_single[['start', 'end']].astype(int)
