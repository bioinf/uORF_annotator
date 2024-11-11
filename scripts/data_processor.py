from logger import Logger
from temporary_file_manager import TemporaryFileManager

from collections import defaultdict
import pandas as pd
import json
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

		return bed_file


class DataProcessor:
	def __init__(self, gtf_file_path, vcf_file_path, bed_file_path, fasta_file_path):
		self.gtf = gtf_file_path
		self.vcf = vcf_file_path
		self.bed = bed_file_path
		self.fasta = fasta_file_path

		self.intersected_bed = ''
		self.bed_4col_info = '|utid|overlapping_type|dominance_type|codon_type'
		self.file_path = ''
		self.bed_4col_info_cols = self.bed_4col_info.split('|')[1:]
		self.source_bed_lines = {}


		self.data = None
		self.first_cds_df = None
		self.uorf_data = None
		self.interorf_single = None
		self.export_exons_file = None
		self.tmp_interorf_single_file = None
		self.source_gtf = None
		self.tmp_io_exon_isec_tab_file = None
		self.interorfs_bed_df = None

	def process_data(self):
		self._initial_process_gtf_file()
		self._logger_0
		self._process_vcf_file()
		self._process_bed_file()
		self._intersect_vcf_and_bed()
		self._load_data(self.intersected_bed)
		self._extract_exon_starts_and_sizes()
		self._extract_predefined_bed_anno()
		self._extract_distancies_from_orf_to_snp()
		self._extract_names()
		self._logger_1()
		self._extract_exons_data()
		self._logger_2()
		self._extract_cds_gtf()
		self._get_first_cds()
		self._extract_uorf_data()
		self._logger_3()
		self._extract_interorf_data()
		self._logger_4()
		self._intersect_interorf_with_exons()
		self._extract_cds_list()
		md5sum_uorf = hashlib.md5(self.uorf_data.to_csv(index=False).encode()).hexdigest()
		md5sum_interorf = hashlib.md5(self.interorf_single.to_csv(index=False).encode()).hexdigest()
		md5sum_interorfs_bed_df = hashlib.md5(self.interorfs_bed_df.to_csv(index=False).encode()).hexdigest()
		print("md5sum:", md5sum_uorf)
		print("md5sum:", md5sum_interorf)
		print("md5sum:", md5sum_interorfs_bed_df)


	def _initial_process_gtf_file(self):
		self.source_gtf = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
		self.source_gtf.columns = ['chr', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'info']
		self.source_gtf['start'] = self.source_gtf['start'].astype(int) - 1
		self.source_gtf['end'] = self.source_gtf['end'].astype(int)
		self.source_gtf = self.source_gtf.loc[self.source_gtf['info'].str.contains('(NM_|ENST)\d+')]
		self.gene_transcript = self.source_gtf['info'].str.extract(r'gene_id \"([^\"]+)', expand=False).to_dict()
		self.source_gtf['transcript'] = self.source_gtf['info'].str.extract(r'transcript_id \"((NM_|ENST)\d+)', expand=False)[0]

	def _process_vcf_file(self):
		self.header_lines = []
		with open(self.vcf, 'r') as f:
			for line in f:
				if line.startswith('##'):
					self.header_lines.append(line.strip())
				elif line.startswith('#CHROM'):
					break
		"\n".join(self.header_lines)

	def _process_bed_file(self):
		orf_count = 1
		with open(self.bed, 'r') as bed_file:
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

	def _intersect_vcf_and_bed(self):
		self.intersected_bed = Bedtools.intersect(self.vcf, self.bed).name

	def _logger_0(self):
		Logger.log_num_rows_and_columns_in_gtf_after_processing(self.source_gtf.shape[0])

	def _logger_1(self):
		Logger.log_num_records_in_table_after_annotation_processing_1(self.data.shape)

	def _logger_2(self):
		Logger.log_num_exons_after_gtf_processing(self.export_exons.shape[0])

	def _logger_3(self):
		Logger.log_num_records_in_table_after_annotation_processing_2(self.uorf_data.shape)

	def _logger_4(self):
		Logger.log_num_records_in_table_after_annotation_processing_3(self.interorf_single.shape)

	def _load_data(self, intersected_bed ):
		self.data = pd.read_table(intersected_bed, header=None)
		self.data = self.data.dropna()
		self.data = self.data.loc[:, [0, 1, 3, 4, 7, 9, 10, 11, 13, 17, 18, 19]]
		self.data = self.data.rename({0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT',
													7: 'INFO', 9: 'orf_start', 10: 'orf_end',
													11: 'bed_anno', 13:'strand', 17:'n_exons', 
													18: 'exons_sizes', 19: 'exons_starts'}, axis=1)

	def _extract_exon_starts_and_sizes(self):
		self.data['exons_starts'] = self.data['exons_starts'].apply(lambda x: list(map(int, x.rstrip(',').split(','))))
		self.data['exons_sizes'] = self.data['exons_sizes'].apply(lambda x: list(map(int, x.rstrip(',').split(','))))

	def _extract_predefined_bed_anno(self):
		self.data = pd.concat([self.data, self.data['bed_anno'].str.split('|', expand=True)], axis=1)
		self.data.columns = list(self.data.columns[:-4]) + self.bed_4col_info_cols

	def _extract_distancies_from_orf_to_snp(self):
		self.data[['POS', 'orf_start', 'orf_end']] = self.data[['POS', 'orf_start', 'orf_end']].astype(int)
		self.data['dist_from_orf_to_snp'] = self.data.apply(lambda x: x['POS'] - x['orf_start'] - 1 if x['strand'] == '+' else x['orf_end'] - x['POS'], axis=1).astype(int)

	def _extract_names(self):
		self.data['name'] = self.data.apply(lambda x: f"{x['#CHROM']}:{x['orf_start']}-{x['orf_end']}({x['strand']})", axis=1)
		self.data['transcript'] = self.data.apply(lambda x: re.findall('[NMENSTX]+_*\d+', x['bed_anno'])[0], axis=1)
		self.data['gene_name'] = self.data.apply(lambda x: self.gene_transcript.get(x['transcript'], ''), axis=1)
		self.data['name_and_trans'] = self.data.apply(lambda x: f"{x['name']}/{x['transcript']}", axis=1)

	def _extract_exons_data(self):
		self.export_exons = self.source_gtf.loc[self.source_gtf['type'] == 'exon']
		self.export_exons = self.export_exons.loc[:, ['chr', 'start', 'end', 'info', 'score', 'strand']]
		self.export_exons = self.export_exons.sort_values(by=['chr', 'start']).drop_duplicates()
		self.export_exons_file = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		self.export_exons.to_csv(self.export_exons_file.name, sep='\t', header=False, index=False)

	def _extract_cds_gtf(self):
		self.cds_df = self.source_gtf.loc[self.source_gtf['type'] == 'CDS']
		self.cds_df = self.cds_df.sort_values(['chr', 'start'])

	def _get_first_cds(self):
		self.first_cds_df = self.cds_df.sort_values(['transcript', 'strand', 'start'])
		self.first_cds_df = self.first_cds_df.groupby('transcript').apply(lambda x: x.iloc[0] if x['strand'].iloc[0] == '+' else x.iloc[-1])
		self.first_cds_df.reset_index(drop=True, inplace=True)
		self.first_cds_df = self.first_cds_df.set_index('transcript')

	def _get_true_cds_start(self, tr_id):
		if tr_id in self.first_cds_df.index:
			return self.first_cds_df.loc[tr_id, 'start'] if self.first_cds_df.loc[tr_id, 'strand'] == '+' else self.first_cds_df.loc[tr_id, 'end']
		return None

	def _extract_uorf_data(self):
		self.uorf_data = self.data.loc[:, ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_starts', 'name', 'name_and_trans']].copy()
		self.uorf_data.columns = ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_norm_starts', 'id', 'name_and_trans']
		self.uorf_data = self.uorf_data.map(lambda x: str(x) if isinstance(x, list) else x).drop_duplicates()
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
			exons_starts = [row['orf_start'] + eStart for eStart in exons_norm_starts]
			return exons_starts, exons_norm_starts
		else:
			exons_starts = [row['orf_start'] + eStart + eSize for eStart, eSize in zip(exons_norm_starts, exons_sizes)]
			exons_norm_starts = [orf_len - eStart - eSize for eStart, eSize in zip(exons_norm_starts, exons_sizes)]
			return exons_starts[::-1], exons_norm_starts[::-1]

	def _extract_interorf_data(self):
		mask = ((self.uorf_data['strand'] == '+') & (self.uorf_data['orf_end'] < self.uorf_data['cds_start'])) | (self.uorf_data['cds_start'] < self.uorf_data['orf_start'])
		interorf_data = self.uorf_data[mask].copy()
		interorf_data.loc[interorf_data['strand'] == '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] == '+', ['orf_end', 'cds_start']].values
		interorf_data.loc[interorf_data['strand'] != '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] != '+', ['cds_start', 'orf_start']].values
		self.interorf_single = interorf_data[['#CHROM', 'start', 'end', 'name_and_trans']].rename(columns={'#CHROM': 'chr'}).dropna()
		self.interorf_single[['start', 'end']] = self.interorf_single[['start', 'end']].astype(int)
		self.tmp_interorf_single_file = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		self.interorf_single.to_csv(self.tmp_interorf_single_file.name, sep='\t', header=False, index=False)


	def _intersect_interorf_with_exons(self):
		self.tmp_io_exon_isec_tab_file = Bedtools.intersect(self.tmp_interorf_single_file.name, self.export_exons_file.name)

		interorfs_bed_dict = {}
		with open(self.tmp_io_exon_isec_tab_file.name, 'r') as isec_handle:
			for line in isec_handle:
				content = line.strip().split('\t')
				uorf_tr_id = content[3].split('/')[1]
				transcript_id_match = re.search(r'transcript_id\s*""([^\.]+)', content[7], re.IGNORECASE).group(1)
				if uorf_tr_id == transcript_id_match:
					rstart = max(int(content[1]), int(content[5]))
					rend = min(int(content[2]), int(content[6]))
					out_line = [content[0], str(rstart), str(rend), content[3], content[8], content[9]]

					key = out_line[3]
					if key not in interorfs_bed_dict:
						interorfs_bed_dict[key] = {'chr': out_line[0], 'strand': out_line[-1], 'exons_sizes': [], 'exons_starts': []}
					interorfs_bed_dict[key]['exons_sizes'].append(int(out_line[2]) - int(out_line[1]))
					interorfs_bed_dict[key]['exons_starts'].append(int(out_line[1]) if out_line[-1] == '+' else int(out_line[2]))

		for k, v in interorfs_bed_dict.items():
			if v['strand'] == '+':
				starts_sorted = sorted(v['exons_starts'])
				v['exons_sizes'] = [x for _, x in sorted(zip(starts_sorted, v['exons_sizes']))]
				v['exons_norm_starts'] = [i-starts_sorted[0] for i in starts_sorted]
			elif v['strand'] == '-':
				size_dict = dict(zip(v['exons_starts'], v['exons_sizes']))
				starts_sorted = sorted(v['exons_starts'], reverse=True)
				v['exons_sizes'] = [size_dict[x] for x in starts_sorted]
				v['exons_norm_starts'] = [abs(i-starts_sorted[0]) for i in starts_sorted]
				v['exons_starts'] = starts_sorted

		self.interorfs_bed_df = pd.DataFrame(interorfs_bed_dict).T
		self.interorfs_bed_df['id'] = self.interorfs_bed_df.index

	def _extract_cds_list(self):
		exons_gtf_list = self.cds_df.values.tolist()
		exons_gtf_dict = defaultdict(lambda: defaultdict(list))
		idx = (0, 3, 4, 6)

		for gtf_line in exons_gtf_list:
			key = re.search('transcript_id \"([^\.]+)', gtf_line[8]).group(1).strip('"')
			values = [gtf_line[i] for i in idx]
			exons_gtf_dict[key]['chr'] = values[0]
			exons_gtf_dict[key]['strand'] = values[-1]
			exons_gtf_dict[key]['exons_sizes'].append(int(values[2])-int(values[1]) + 1)
			exons_gtf_dict[key]['exons_starts'].append(int(values[1]) - 1 if values[-1] == '+' else int(values[2]))

		for k in exons_gtf_dict.keys():
			if exons_gtf_dict[k]['strand'] == '+':
				starts_sorted = sorted(exons_gtf_dict[k]['exons_starts'])
				exons_gtf_dict[k]['exons_sizes'] = [x for _, x in sorted(zip(starts_sorted, exons_gtf_dict[k]['exons_sizes']))]
				exons_gtf_dict[k]['id'] = f"{exons_gtf_dict[k]['chr']}:{starts_sorted[0]}(+)"
				exons_gtf_dict[k]['exons_norm_starts'] = [i-starts_sorted[0] for i in starts_sorted]
			elif exons_gtf_dict[k]['strand'] == '-':
				size_dict = {k: v for k, v in zip(exons_gtf_dict[k]['exons_starts'], exons_gtf_dict[k]['exons_sizes'])}
				starts_sorted = sorted(exons_gtf_dict[k]['exons_starts'], reverse=True)
				exons_gtf_dict[k]['exons_starts'] = starts_sorted
				exons_gtf_dict[k]['exons_sizes'] = [size_dict[x] for x in starts_sorted]
				exons_gtf_dict[k]['id'] = f"{exons_gtf_dict[k]['chr']}:{starts_sorted[0]}(-)"
				exons_gtf_dict[k]['exons_norm_starts'] = [abs(i-starts_sorted[0]) for i in starts_sorted]
		# print(exons_gtf_dict)
		json.dumps(exons_gtf_dict)
		exons_gtf_dict = hashlib.md5(json.dumps(exons_gtf_dict).encode()).hexdigest()
		print("md5sum:", exons_gtf_dict)
