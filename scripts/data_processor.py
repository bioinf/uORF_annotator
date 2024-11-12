from typing import Dict, List, Optional, Tuple

from logger import Logger
from temporary_file_manager import TemporaryFileManager

from collections import defaultdict
from pybedtools import BedTool
from Bio import SeqIO
import pandas as pd
import json
import re

import hashlib

class Bedtools:
	@staticmethod
	def intersect(file1, file2):
		bed_file1 = BedTool(file1)
		bed_file2 = BedTool(file2)
		intersected_bed = bed_file1.intersect(bed_file2, wo=True)
		num_lines = len(intersected_bed)
		Logger.log_num_variants_in_intersection(num_lines)
		intersection_df = intersected_bed.to_dataframe(header=None)
		intersection_df.columns = range(len(intersection_df.columns))
		return intersection_df
	
	@staticmethod
	def getfasta(fasta_file, bed_file, strand=True):
		fasta = BedTool(fasta_file)
		bed = BedTool(bed_file)
		getfasta = bed.sequence(fi=fasta, s=strand)
		return getfasta.seqfn


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
		self.raw_bed_lines = {}


		self.uorf_variation_data = None
		self.first_cds_df = None
		self.uorf_data = None
		self.interorf_single = None
		self.exons_data_file = None
		self.tmp_interorf_single_file = None
		self.gtf_raw = None
		self.tmp_io_exon_isec_tab_file = None
		self.interorfs_bed_df = None

	def process_data(self):
		self.gtf_raw, self.gene_transcript = self._initial_process_gtf_file(self.gtf)
		self.vcf_header = self._process_vcf_file(self.vcf)
		self.raw_bed_lines = self._process_bed_file(self.bed)
		self.intersected_bed = self._intersect_vcf_and_bed(self.vcf, self.bed)

		self.uorf_variation_data = self._load_data(self.intersected_bed)
		self.uorf_variation_data = self._extract_uorf_exon_bounds(self.uorf_variation_data)
		self.uorf_variation_data = self._extract_uorf_annotations(self.uorf_variation_data, self.bed_4col_info_cols)
		self.uorf_variation_data = self._calculate_uorf_snp_distances(self.uorf_variation_data)
		self.uorf_variation_data = self._extract_gene_info(self.uorf_variation_data)

		self.exons_data, self.exons_data_file = self._extract_exons_data(self.gtf_raw)
		self.cds_df = self._extract_cds_gtf(self.gtf_raw)
		self.first_cds_df = self._get_first_cds(self.cds_df)
		self.uorf_data = self._extract_uorf_data(self.uorf_variation_data, self.first_cds_df)
		self.interorf_single, self.tmp_interorf_single_file = self._extract_interorf_data(self.uorf_data)
		self.interorfs_bed_df, self.tmp_interorfs_bed_sorted = self._intersect_interorf_with_exons(self.tmp_interorf_single_file, self.exons_data_file)
		self.exons_gtf_dict = self._extract_cds_list(self.cds_df)
		self.uorf_dict = self._get_uorf_dict(self.fasta, self.bed)
		self.interorfs_dict = self._get_interorfs_fasta(self.fasta, self.tmp_interorfs_bed_sorted)
		self.cds_dict = self._extract_cds_regions(self.fasta, self.cds_df)
	
		self._logger_0
		self._logger_1()
		self._logger_2()
		self._logger_3()
		self._logger_4()
		md5sum_uorf = hashlib.md5(self.uorf_data.to_csv(index=False).encode()).hexdigest()
		md5sum_interorf = hashlib.md5(self.interorf_single.to_csv(index=False).encode()).hexdigest()
		md5sum_interorfs_bed_df = hashlib.md5(self.interorfs_bed_df.to_csv(index=False).encode()).hexdigest()
		exons_gtf_dict = hashlib.md5(json.dumps(self.exons_gtf_dict).encode()).hexdigest()
		print("md5sum:", md5sum_uorf)
		print("md5sum:", md5sum_interorf)
		print("md5sum:", md5sum_interorfs_bed_df)
		print("md5sum:", exons_gtf_dict)
		md5sum = hashlib.md5(str(self.interorfs_dict).encode()).hexdigest()
		print(md5sum)
		md5sum = hashlib.md5(str(self.cds_dict).encode()).hexdigest()
		print(md5sum)

	def _get_uorf_dict(self, fasta, bed, strand=True):
		getfasta = Bedtools.getfasta(fasta, bed, strand)
		uorf_dict = SeqIO.to_dict(SeqIO.parse(getfasta, "fasta"))
		return uorf_dict

	def _initial_process_gtf_file(self, gtf_path: str) -> Tuple[pd.DataFrame, Dict[str, str]]:
		gtf_raw = pd.read_csv(gtf_path, sep='\t', comment='#', header=None)
		gtf_raw.columns = ['chr', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'info']
		gtf_raw['start'] = gtf_raw['start'].astype(int) - 1
		gtf_raw['end'] = gtf_raw['end'].astype(int)
		gtf_raw = gtf_raw.loc[gtf_raw['info'].str.extract('(NM_|ENST)\\d+').any(axis=1)]
		gene_transcript = gtf_raw['info'].str.extract(r'gene_id \"([^\"]+)', expand=False).to_dict()
		gtf_raw['transcript'] = gtf_raw['info'].str.extract(r'transcript_id \"((NM_|ENST)\d+)', expand=False)[0]
		return gtf_raw, gene_transcript

	def _process_vcf_file(self, vcf_path: str) -> List[str]:
		vcf_header = []
		with open(vcf_path, 'r') as f:
			for line in f:
				if line.startswith('##'):
					vcf_header.append(line.strip())
				elif line.startswith('#CHROM'):
					break
		
		return vcf_header

	def _process_bed_file(self, bed_path: str) -> Dict[str, List[str]]:
		raw_bed_lines = {}
		orf_count = 1
		
		with open(bed_path, 'r') as bed_file:
			for line in bed_file:
				if line.startswith('#'):
					continue
				columns = line.strip().split('\t')
				trans_id = self._extract_trans_id(columns[3])
				uo_name = self._create_uo_name(columns, trans_id)
				columns[3] = self._update_custom_annotation(orf_count, columns[3])
				orf_count += 1
				raw_bed_lines[uo_name] = columns.copy()
				
		return raw_bed_lines

	def _extract_trans_id(self, column_3):
		return re.findall('[NMENST]+_*\d+', column_3)[0]

	def _create_uo_name(self, columns, trans_id):
		return f"{columns[0]}:{columns[1]}-{columns[2]}({columns[5]})/{trans_id}"

	def _update_custom_annotation(self, orf_count, column_3):
		_, overlap_type, _, codon_type = column_3.split('|')
		return f"UORF{orf_count}|{overlap_type}|{codon_type}"

	def get_raw_bed_lines(self):
		return self.raw_bed_lines

	def _intersect_vcf_and_bed(self, vcf_path: str, bed_path: str) -> pd.DataFrame:
		intersected_bed = Bedtools.intersect(vcf_path, bed_path)
		# Logger.log_num_variants_in_intersection(
		# 	int(subprocess.check_output(f"wc -l {intersected_bed.name}", shell=True).split()[0])
		# )
		return intersected_bed

	def _logger_0(self):
		Logger.log_num_rows_and_columns_in_gtf_after_processing(self.gtf_raw.shape[0])

	def _logger_1(self):
		Logger.log_num_records_in_table_after_annotation_processing_1(self.uorf_variation_data.shape)

	def _logger_2(self):
		Logger.log_num_exons_after_gtf_processing(self.exons_data.shape[0])

	def _logger_3(self):
		Logger.log_num_records_in_table_after_annotation_processing_2(self.uorf_data.shape)

	def _logger_4(self):
		Logger.log_num_records_in_table_after_annotation_processing_3(self.interorf_single.shape)

	def _load_data(self, uorf_variants: pd.DataFrame) -> pd.DataFrame:
		uorf_variation_data = uorf_variants.dropna()
	
		# Select and rename relevant columns
		columns_to_keep = [0, 1, 3, 4, 7, 9, 10, 11, 13, 17, 18, 19]
		column_names = {
			0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT',
			7: 'INFO', 9: 'orf_start', 10: 'orf_end',
			11: 'bed_anno', 13:'strand', 17:'n_exons', 
			18: 'exons_sizes', 19: 'exons_starts'
		}
	
		uorf_variation_data = uorf_variation_data.loc[:, columns_to_keep]
		uorf_variation_data = uorf_variation_data.rename(column_names, axis=1)
	
		return uorf_variation_data

	def _extract_uorf_exon_bounds(self, uorf_variation_data: pd.DataFrame) -> pd.DataFrame:
		uorf_exon_bounds = uorf_variation_data.copy()
		uorf_exon_bounds['exons_starts'] = uorf_exon_bounds['exons_starts'].apply(
			lambda x: list(map(int, x.rstrip(',').split(',')))
		)
		uorf_exon_bounds['exons_sizes'] = uorf_exon_bounds['exons_sizes'].apply(
			lambda x: list(map(int, x.rstrip(',').split(',')))
		)
		return uorf_exon_bounds

	def _extract_uorf_annotations(self, uorf_variation_data: pd.DataFrame, bed_annotation_columns: List[str]) -> pd.DataFrame:
		uorf_anno = uorf_variation_data.copy()
		bed_anno_split = uorf_anno['bed_anno'].str.split('|', expand=True)
		uorf_info_data = pd.concat([uorf_anno, bed_anno_split], axis=1)
		uorf_info_data.columns = list(uorf_info_data.columns[:-4]) + bed_annotation_columns
		return uorf_info_data

	def _calculate_uorf_snp_distances(self, uorf_variation_data: pd.DataFrame) -> pd.DataFrame:
		uorf_variation_data[['POS', 'orf_start', 'orf_end']] = uorf_variation_data[['POS', 'orf_start', 'orf_end']].astype(int)
		uorf_variation_data['dist_from_orf_to_snp'] = uorf_variation_data.apply(lambda x: x['POS'] - x['orf_start'] - 1 if x['strand'] == '+' else x['orf_end'] - x['POS'], axis=1).astype(int)
		return uorf_variation_data

	def _extract_gene_info(self, uorf_variation_data: pd.DataFrame) -> pd.DataFrame:
		uorf_variation_data['name'] = uorf_variation_data.apply(lambda x: f"{x['#CHROM']}:{x['orf_start']}-{x['orf_end']}({x['strand']})", axis=1)
		uorf_variation_data['transcript'] = uorf_variation_data.apply(lambda x: re.findall('[NMENSTX]+_*\d+', x['bed_anno'])[0], axis=1)
		uorf_variation_data['gene_name'] = uorf_variation_data.apply(lambda x: self.gene_transcript.get(x['transcript'], ''), axis=1)
		uorf_variation_data['name_and_trans'] = uorf_variation_data.apply(lambda x: f"{x['name']}/{x['transcript']}", axis=1)
		return uorf_variation_data

	def _extract_exons_data(self, gtf_raw: pd.DataFrame) -> pd.DataFrame:
		exons_data = gtf_raw.loc[gtf_raw['type'] == 'exon']
		exons_data = exons_data.loc[:, ['chr', 'start', 'end', 'info', 'score', 'strand']]
		exons_data = exons_data.sort_values(by=['chr', 'start']).drop_duplicates()
		exons_data_file = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		exons_data.to_csv(exons_data_file.name, sep='\t', header=False, index=False)
		return exons_data, exons_data_file

	def _extract_cds_gtf(self, gtf_raw):
		cds_df = gtf_raw.loc[gtf_raw['type'] == 'CDS']
		cds_df = cds_df.sort_values(by=['chr', 'start']).drop_duplicates()
		return cds_df

	def _get_first_cds(self, cds_df):
		first_cds_df = cds_df.sort_values(['transcript', 'strand', 'start'])
		first_cds_df = first_cds_df.groupby('transcript').apply(lambda x: x.iloc[0] if x['strand'].iloc[0] == '+' else x.iloc[-1])
		first_cds_df.reset_index(drop=True, inplace=True)
		first_cds_df = first_cds_df.set_index('transcript')
		return first_cds_df

	def _get_true_cds_start(self, tr_id, first_cds_df):
		if tr_id in first_cds_df.index:
			return first_cds_df.loc[tr_id, 'start'] if first_cds_df.loc[tr_id, 'strand'] == '+' else first_cds_df.loc[tr_id, 'end']
		return None

	def _extract_uorf_data(self, uorf_variation_data: pd.DataFrame, first_cds_df: pd.DataFrame) -> pd.DataFrame:
		uorf_data = uorf_variation_data.loc[:, ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_starts', 'name', 'name_and_trans']].copy()
		uorf_data.columns = ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_norm_starts', 'id', 'name_and_trans']
		uorf_data = uorf_data.map(lambda x: str(x) if isinstance(x, list) else x).drop_duplicates()
		uorf_data['exons_norm_starts'] = uorf_data['exons_norm_starts'].apply(eval)
		uorf_data['exons_sizes'] = uorf_data['exons_sizes'].apply(eval)
		uorf_data[['exons_starts', 'exons_norm_starts']] = uorf_data.apply(self._process_exons, axis=1, result_type='expand')
		uorf_data['exons_sizes'] = uorf_data.apply(lambda row: row['exons_sizes'] if row['strand'] == '+' else row['exons_sizes'][::-1], axis=1)
		uorf_data.index = uorf_data['id']
		uorf_data['cds_start'] = uorf_data['transcript'].apply(lambda x: self._get_true_cds_start(x, first_cds_df))
		return uorf_data

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

	def _extract_interorf_data(self, uorf_data: pd.DataFrame):
		mask = ((uorf_data['strand'] == '+') & (uorf_data['orf_end'] < uorf_data['cds_start'])) | (uorf_data['cds_start'] < uorf_data['orf_start'])
		interorf_data = uorf_data[mask].copy()
		interorf_data.loc[interorf_data['strand'] == '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] == '+', ['orf_end', 'cds_start']].values
		interorf_data.loc[interorf_data['strand'] != '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] != '+', ['cds_start', 'orf_start']].values
		interorf_single = interorf_data[['#CHROM', 'start', 'end', 'name_and_trans']].rename(columns={'#CHROM': 'chr'}).dropna()
		interorf_single[['start', 'end']] = interorf_single[['start', 'end']].astype(int)
		tmp_interorf_single_file = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		interorf_single.to_csv(tmp_interorf_single_file.name, sep='\t', header=False, index=False)
		return interorf_single, tmp_interorf_single_file


	def _intersect_interorf_with_exons(self, tmp_interorf_single_file, exons_data_file):
		tmp_io_exon_isec_tab_file = Bedtools.intersect(tmp_interorf_single_file.name, exons_data_file.name)
		tmp_io_exon_isec_tab_file.columns = ['chr', 'start', 'end', 'id', 'chr_exon', 'start_exon', 'end_exon', 'gene_info', 'score', 'strand', 'length']
		interorfs_bed_dict = {}

		sorted_lines = []
		for _, row in tmp_io_exon_isec_tab_file.iterrows():
			uorf_tr_id = row['id'].split('/')[1]
			transcript_id_match = re.search(r'transcript_id\s*"([^\.]+)', row['gene_info'], re.IGNORECASE).group(1)
			if uorf_tr_id == transcript_id_match:
				rstart = max(int(row['start']), int(row['start_exon']))
				rend = min(int(row['end']), int(row['end_exon']))
				out_line = [row['chr'], str(rstart), str(rend), row['id'], row['score'], row['strand']]

				key = out_line[3]
				if key not in interorfs_bed_dict:
					interorfs_bed_dict[key] = {'chr': out_line[0], 'strand': out_line[-1], 'exons_sizes': [], 'exons_starts': []}
				interorfs_bed_dict[key]['exons_sizes'].append(int(out_line[2]) - int(out_line[1]))
				interorfs_bed_dict[key]['exons_starts'].append(int(out_line[1]) if out_line[-1] == '+' else int(out_line[2]))
				sorted_lines.append(out_line)

		tmp_interorfs_bed_sorted = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		with open(tmp_interorfs_bed_sorted.name, 'w') as bed_file:
			for line in sorted(sorted_lines, key=lambda x: (x[0], int(x[1]))):
				bed_file.write('\t'.join(line) + '\n')

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

		interorfs_bed_df = pd.DataFrame(interorfs_bed_dict).T
		interorfs_bed_df['id'] = interorfs_bed_df.index

		return interorfs_bed_df, tmp_interorfs_bed_sorted

	def _extract_cds_list(self, cds_df):
		exons_gtf_list = cds_df.values.tolist()
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

		return exons_gtf_dict

	def _get_interorfs_fasta(self, fasta, tmp_interorfs_bed_sorted, strand=True):
		tmp_interorfs_fasta = Bedtools.getfasta(fasta, tmp_interorfs_bed_sorted.name, strand)

		tmp_interorfs_full_seq = TemporaryFileManager.create('.bed', tmp_dir='/tmp')

		seen_parts = set()
		interorfs_full_seq = defaultdict(str)

		with open(tmp_interorfs_fasta, "r") as split_fasta:
			current_name = None
			for line in split_fasta:
				if line.startswith('>'):
					record_tmp_name = line[1:].strip()
					uorf_tmp_name = record_tmp_name.split('::')[0]

					if record_tmp_name in seen_parts:
						current_name = None
					else:
						seen_parts.add(record_tmp_name)
						current_name = uorf_tmp_name
				elif current_name:
					if '(-)' in current_name:
						interorfs_full_seq[current_name] = line.strip() + interorfs_full_seq[current_name]
					else:
						interorfs_full_seq[current_name] += line.strip()

		with open(tmp_interorfs_full_seq.name, 'w') as full_fasta:
			for uorf_name, seq in interorfs_full_seq.items():
				full_fasta.write(f">{uorf_name}\n{seq}\n")

		interorfs_dict = SeqIO.to_dict(SeqIO.parse(tmp_interorfs_full_seq.name, "fasta"))

		return interorfs_dict

	def _extract_cds_regions(self, fasta, cds_df):
		cds_bed_df = cds_df[['chr', 'start', 'end', 'strand']]
		cds_bed_df.loc[:, 'start'] = cds_bed_df['start'].astype(int, copy=False)
		cds_bed_df.loc[:, 'end'] = cds_bed_df['end'].astype(int, copy=False)
		cds_bed_df.loc[:, 'start'] -= 1
		cds_bed_df = cds_bed_df.sort_values(by=['chr', 'start'])
		cds_bed_df = cds_bed_df.drop_duplicates()

		cds_bed = TemporaryFileManager.create('.bed', tmp_dir='/tmp')

		cds_bed_df.to_csv(cds_bed.name, sep='\t', index=False, header=False)

		getfasta = Bedtools.getfasta(fasta, cds_bed.name, strand=True)

		cds_dict = SeqIO.to_dict(SeqIO.parse(getfasta, "fasta"))

		return cds_dict

	def _create_complete_cds(self, chrom, strand, eStarts, eSizes):
		if chrom in ['chrM', 'chrMT', 'M', 'MT']:
			return None
		out_seq = ''
		for i, j in zip(eStarts, eSizes):
			if strand == '+':
				out_seq += self.cds_dict[f'{chrom}:{i}-{i+j}({strand})'].seq
			else:
				out_seq += self.cds_dict[f'{chrom}:{i-j}-{i}({strand})'].seq
		return out_seq

	def _process_exons_gtf_dict(self, exons_gtf_dict):
		for k in exons_gtf_dict.keys():
			echrom = exons_gtf_dict[k]['chr']
			estrand = exons_gtf_dict[k]['strand']
			estarts = exons_gtf_dict[k]['exons_starts']
			esizes = exons_gtf_dict[k]['exons_sizes']
			exons_gtf_dict[k]['seq'] = str(self._create_complete_cds(echrom, estrand, estarts, esizes))
		exons_gtf_df = pd.DataFrame(exons_gtf_dict).T
		return exons_gtf_df

