from logger import Logger
from temporary_file_manager import TemporaryFileManager

import re
import pandas as pd

class BEDFile:

	def __init__(self, file_path):
		self.file_path = file_path
		self.bed_4col_info = '|utid|overlapping_type|dominance_type|codon_type'
		self.bed_4col_info_cols = self.bed_4col_info.split('|')[1:]
		self.source_bed_lines = {}

	def process_bed_file(self):
		self._parse_file()
		self._logger()

	def _parse_file(self):
		orf_count = 1
		with open(self.file_path, 'r') as bed_file:
			for line in bed_file:
				if line.startswith('#'):
					continue
				columns = line.strip().split('\t')
				trans_id = self._extract_trans_id(columns[3])
				uo_name = self._create_uo_name(columns, trans_id)
				columns[3] = self._update_custom_annotation(orf_count, columns[3])
				orf_count += 1
				self.source_bed_lines[uo_name] = columns.copy()

	def _extract_trans_id(self, column_3):
		return re.findall('[NMENST]+_*\d+', column_3)[0]

	def _create_uo_name(self, columns, trans_id):
		return f"{columns[0]}:{columns[1]}-{columns[2]}({columns[5]})/{trans_id}"

	def _update_custom_annotation(self, orf_count, column_3):
		_, overlap_type, _, codon_type = column_3.split('|')
		return f"UORF{orf_count}|{overlap_type}|{codon_type}"

	def get_source_bed_lines(self):
		return self.source_bed_lines

	def _logger(self):
		Logger.log_processed_bed_records(len(self.source_bed_lines))
