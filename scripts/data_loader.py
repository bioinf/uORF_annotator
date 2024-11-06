from logger import Logger
from temporary_file_manager import TemporaryFileManager

import re
import pandas as pd


class GTFFile:
	def __init__(self, file_path):
		self.file_path = file_path
		self.source_gtf = None
		self.gene_transcript = None

	def process_gtf_file(self):
		self._load_gtf_data()
		self._extract_gene_transcript()
		self._extract_bed_anno()
		self._logger()

	def _load_gtf_data(self):
		self.source_gtf = pd.read_csv(self.file_path, sep='\t', comment='#', header=None)
		self.source_gtf.columns = ['chr', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'info']
		self.source_gtf['start'] = self.source_gtf['start'].astype(int) - 1
		self.source_gtf['end'] = self.source_gtf['end'].astype(int)
		self.source_gtf = self.source_gtf.loc[self.source_gtf['info'].str.contains('(NM_|ENST)\d+')]

	def _extract_gene_transcript(self):
		self.gene_transcript = self.source_gtf['info'].str.extract(r'gene_id \"([^\"]+)', expand=False).to_dict()

	def _extract_bed_anno(self):
		self.source_gtf['transcript'] = self.source_gtf['info'].str.extract(r'transcript_id \"((NM_|ENST)\d+)', expand=False)[0]

	def _logger(self):
		Logger.log_num_rows_and_columns_in_gtf_after_processing(self.source_gtf.shape[0])

class VCFFile:
	def __init__(self, file_path):
		self.file_path = file_path
		self.header_lines = []

	def process_vcf_file(self):
		self._read_header()

	def _read_header(self):
		with open(self.file_path, 'r') as f:
			for line in f:
				if line.startswith('##'):
					self.header_lines.append(line.strip())
				elif line.startswith('#CHROM'):
					break

		return "\n".join(self.header_lines)

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
