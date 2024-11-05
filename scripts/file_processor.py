from logger import Logger
from temporary_file_manager import TemporaryFileManager

import re
import gzip
import pandas as pd


class GTFFile:
	def __init__(self, file_path):
		self.file_path = file_path

	def get_gene_transcript_records(self):
		gene_transcript = {}
		with gzip.open(self.file_path, 'rt') as gtf_handle:
			for line in gtf_handle:
				if 'gene_id' not in line and 'transcript_id' not in line:
					continue
				gene_name = re.findall('gene_id \"([^\"]+)', line)[0]
				transcript_id = re.findall('transcript_id \"([^\"^\.]+)', line)[0]
				gene_transcript[transcript_id] = gene_name
		Logger.log_gene_transcript_records(len(gene_transcript))
		return gene_transcript

	def process_gtf_file(self):

		source_gtf = pd.read_csv(self.file_path, sep='\t', comment='#', header=None)
		source_gtf.columns = ['chr', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'info']
		source_gtf['start'] = [int(x) - 1 for x in source_gtf['start']]
		source_gtf['end'] = [int(x) for x in source_gtf['end']]
		ok_rows = [bool(re.search('(NM_|ENST)\d+', x)) for x in source_gtf['info']]
		source_gtf = source_gtf.loc[ok_rows]
		source_gtf['transcript'] = [re.findall('transcript_id \"((NM_|ENST)\d+)', x)[0][0] for x in source_gtf['info']]

		exons_gtf = source_gtf.loc[source_gtf['type'] == 'exon']
		export_exons = exons_gtf.loc[:, ['chr', 'start', 'end', 'transcript', 'score', 'strand']]
		export_exons = export_exons.sort_values(by=['chr', 'transcript'])
		export_exons = export_exons.drop_duplicates()

		Logger.log_num_exons_after_gtf_processing(export_exons.shape[0])


		tmp_exons_bed = TemporaryFileManager.create('.bed')
		export_exons.to_csv(tmp_exons_bed.name, sep='\t', header=False, index=False)

		tmp_exons_bed.register_at_exit()

		return tmp_exons_bed

class VCFFile:
	def __init__(self, file_path):
		self.file_path = file_path

	def read_header(self):
		header_lines = []
		with open(self.file_path, 'r') as f:
			for line in f:
				if line.startswith('##'):
					header_lines.append(line.strip())
				elif line.startswith('#CHROM'):
					break

		return "\n".join(header_lines)

class BEDFile:
	bed_4col_info = '|utid|overlapping_type|dominance_type|codon_type'
	bed_4col_info_cols = bed_4col_info.split('|')[1:]

	def __init__(self, file_path):
		self.file_path = file_path

	def process_bed_file(self):
		orf_count = 1
		with open(self.file_path, 'r') as bed_file:
			source_bed_lines = {}
			for line in bed_file:
				if line.startswith('#'):
					continue
				columns = line.strip().split('\t')
				trans_id = self._extract_trans_id(columns[3])
				uo_name = self._create_uo_name(columns, trans_id)
				columns[3] = self._update_custom_annotation(orf_count, columns[3])
				orf_count += 1
				source_bed_lines[uo_name] = columns.copy()
		Logger.log_processed_bed_records(len(source_bed_lines))
		return source_bed_lines

	def _extract_trans_id(self, column_3):
		return re.findall('[NMENST]+_*\d+', column_3)[0]

	def _create_uo_name(self, columns, trans_id):
		return f"{columns[0]}:{columns[1]}-{columns[2]}({columns[5]})/{trans_id}"

	def _update_custom_annotation(self, orf_count, column_3):
		_, overlap_type, _, codon_type = column_3.split('|')
		return f"UORF{orf_count}|{overlap_type}|{codon_type}"

	def get_source_bed_lines(self):
		return self.source_bed_lines
