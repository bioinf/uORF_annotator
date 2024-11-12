from logger import Logger
from temporary_file_manager import TemporaryFileManager
from bedtools_wrapper import Bedtools

from typing import Dict, List, Tuple, Optional

import re
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class RawDataProcessor:
	def __init__(self):
		self.vcf_file_path: str = ''
		self.bed_file_path: str = ''
		self.gtf_file_path: str = ''
		self.fasta_file_path: str = ''

		self.vcf_header: List[str] = None
		self.raw_bed_lines: Dict[str, List[str]] = None
		self.gene_transcript: Dict[str, str] = None
		self.gtf_raw: pd.DataFrame = None

	def process_data(self, vcf_file_path: str, bed_file_path: str, gtf_file_path: str, fasta_file_path: str) -> None:

		self.vcf_file_path = vcf_file_path
		self.bed_file_path = bed_file_path
		self.gtf_file_path = gtf_file_path
		self.fasta_file_path = fasta_file_path

		self.process_vcf_file()
		self.process_bed_file()
		self.process_gtf_file()
		self.process_fasta_file()

	def process_vcf_file(self) -> None:
		self.vcf_header = self._process_vcf_header(self.vcf_file_path)
	
	def process_bed_file(self) -> None:
		self.raw_bed_lines = self._process_bed_file(self.bed_file_path)

	def process_gtf_file(self) -> None:
		self.gtf_raw, self.gene_transcript = self._initial_process_gtf_file(self.gtf_file_path)
	
	def process_fasta_file(self) -> None:
		self.fasta = self.fasta_file_path

	def _process_vcf_header(self, vcf_path: str) -> List[str]:
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
				trans_id = self._get_trans_id(columns[3])
				uo_name = self._create_uo_name(columns, trans_id)
				columns[3] = self._update_custom_annotation(orf_count, columns[3])
				orf_count += 1
				raw_bed_lines[uo_name] = columns.copy()
		Logger.log_num_bed_lines(len(raw_bed_lines))
		return raw_bed_lines

	def _get_trans_id(self, column_3: str) -> str:
		return re.findall('[NMENST]+_*\d+', column_3)[0]

	def _create_uo_name(self, columns: List[str], trans_id: str) -> str:
		return f"{columns[0]}:{columns[1]}-{columns[2]}({columns[5]})/{trans_id}"

	def _update_custom_annotation(self, orf_count: int, column_3: str) -> str:
		_, overlap_type, _, codon_type = column_3.split('|')
		return f"UORF{orf_count}|{overlap_type}|{codon_type}"

	def _initial_process_gtf_file(self, gtf_path: str) -> Tuple[pd.DataFrame, Dict[str, str]]:
		gtf_raw = pd.read_csv(gtf_path, sep='\t', comment='#', header=None)
		gtf_raw.columns = ['chr', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'info']
		gtf_raw['start'] = gtf_raw['start'].astype(int) - 1
		gtf_raw['end'] = gtf_raw['end'].astype(int)
		gtf_raw = gtf_raw.loc[gtf_raw['info'].str.extract('(NM_|ENST)\\d+').any(axis=1)]
		gene_transcript = gtf_raw['info'].str.extract(r'gene_id \"([^\"]+)', expand=False).to_dict()
		gtf_raw['transcript'] = gtf_raw['info'].str.extract(r'transcript_id \"((NM_|ENST)\d+)', expand=False)[0]
		Logger.log_num_gene_transcript_records(len(gene_transcript))
		return gtf_raw, gene_transcript


class uORFDataProcessor:
	def __init__(self):
		self.uorf_dict: Dict[str, SeqRecord] = None
		self.intersected_bed: pd.DataFrame = None
		self.uorf_variation_data: pd.DataFrame = None

	def process_data(self, bed_4col_info_cols: List[str], raw_data_processor: RawDataProcessor) -> None:

		self.bed_4col_info_cols = bed_4col_info_cols
		self.vcf_file_path = raw_data_processor.vcf_file_path
		self.bed_file_path = raw_data_processor.bed_file_path
		self.fasta_file_path = raw_data_processor.fasta_file_path
		self.gene_transcript = raw_data_processor.gene_transcript

		self.process_uorf_data()


	def process_uorf_data(self) -> None:
		self.uorf_dict = self._get_uorf_dict(self.fasta_file_path, self.bed_file_path)
		self.intersected_bed = self._get_variation_data_from_known_uorfs(self.vcf_file_path, self.bed_file_path)
		self.uorf_variation_data = self._prepare_uorf_variation_data(self.intersected_bed)
		self.uorf_variation_data = self._get_uorf_exon_bounds(self.uorf_variation_data)
		self.uorf_variation_data = self._get_uorf_annotations(self.uorf_variation_data, self.bed_4col_info_cols)
		self.uorf_variation_data = self._calculate_uorf_snp_distances(self.uorf_variation_data)
		self.uorf_variation_data = self._get_gene_info(self.uorf_variation_data, self.gene_transcript)

	def _get_uorf_dict(self, fasta: str, bed: str, strand: bool = True) -> Dict[str, SeqRecord]:
		getfasta = Bedtools.getfasta(fasta, bed, strand)
		uorf_dict = SeqIO.to_dict(SeqIO.parse(getfasta, "fasta"))
		Logger.log_num_uorf_dict_records(len(uorf_dict))
		return uorf_dict

	def _get_variation_data_from_known_uorfs(self, vcf_path: str, bed_path: str) -> pd.DataFrame:
		intersected_bed = Bedtools.intersect(vcf_path, bed_path)
		Logger.log_num_intersected_bed_records(len(intersected_bed))
		return intersected_bed

	def _prepare_uorf_variation_data(self, uorf_variants: pd.DataFrame) -> pd.DataFrame:
		uorf_variation_data = uorf_variants.dropna()
		columns_to_keep = [0, 1, 3, 4, 7, 9, 10, 11, 13, 17, 18, 19]
		column_names = {
			0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT',
			7: 'INFO', 9: 'orf_start', 10: 'orf_end',
			11: 'bed_anno', 13:'strand', 17:'n_exons', 
			18: 'exons_sizes', 19: 'exons_starts'
		}
		uorf_variation_data = uorf_variation_data.loc[:, columns_to_keep]
		uorf_variation_data = uorf_variation_data.rename(column_names, axis=1)
		Logger.log_num_uorf_variation_data_records(uorf_variation_data.shape)
		return uorf_variation_data

	def _get_uorf_exon_bounds(self, uorf_variation_data: pd.DataFrame) -> pd.DataFrame:
		uorf_exon_bounds = uorf_variation_data.copy()
		uorf_exon_bounds['exons_starts'] = uorf_exon_bounds['exons_starts'].apply(
			lambda x: list(map(int, x.rstrip(',').split(',')))
		)
		uorf_exon_bounds['exons_sizes'] = uorf_exon_bounds['exons_sizes'].apply(
			lambda x: list(map(int, x.rstrip(',').split(',')))
		)
		Logger.log_num_uorf_exon_bounds_records(uorf_variation_data.shape)
		return uorf_exon_bounds

	def _get_uorf_annotations(self, uorf_variation_data: pd.DataFrame, bed_annotation_columns: List[str]) -> pd.DataFrame:
		uorf_anno = uorf_variation_data.copy()
		bed_anno_split = uorf_anno['bed_anno'].str.split('|', expand=True)
		uorf_info_data = pd.concat([uorf_anno, bed_anno_split], axis=1)
		uorf_info_data.columns = list(uorf_info_data.columns[:-4]) + bed_annotation_columns
		Logger.log_num_uorf_annotations_records(uorf_variation_data.shape)
		return uorf_info_data

	def _calculate_uorf_snp_distances(self, uorf_variation_data: pd.DataFrame) -> pd.DataFrame:
		uorf_variation_data[['POS', 'orf_start', 'orf_end']] = uorf_variation_data[['POS', 'orf_start', 'orf_end']].astype(int)
		uorf_variation_data['dist_from_orf_to_snp'] = uorf_variation_data.apply(lambda x: x['POS'] - x['orf_start'] - 1 if x['strand'] == '+' else x['orf_end'] - x['POS'], axis=1).astype(int)
		Logger.log_num_uorf_snp_distances_records(uorf_variation_data.shape)
		return uorf_variation_data
	
	def _get_gene_info(self, uorf_variation_data: pd.DataFrame, gene_transcript: Dict[str, str] = {}) -> pd.DataFrame:
		uorf_variation_data['name'] = uorf_variation_data.apply(lambda x: f"{x['#CHROM']}:{x['orf_start']}-{x['orf_end']}({x['strand']})", axis=1)
		uorf_variation_data['transcript'] = uorf_variation_data.apply(lambda x: re.findall('[NMENSTX]+_*\d+', x['bed_anno'])[0], axis=1)
		uorf_variation_data['gene_name'] = uorf_variation_data.apply(lambda x: gene_transcript.get(x['transcript'], ''), axis=1)
		uorf_variation_data['name_and_trans'] = uorf_variation_data.apply(lambda x: f"{x['name']}/{x['transcript']}", axis=1)
		Logger.log_num_gene_info_records(uorf_variation_data.shape)
		return uorf_variation_data
	
class InterORFDataProcessor:
	def __init__(self):
		self.exons_data: pd.DataFrame = None
		self.exons_data_file: str = None
		self.cds_df: pd.DataFrame = None
		self.exons_gtf_dict: Dict[str, Dict[str, List[int]]] = None
		self.first_cds_df: pd.DataFrame = None
		self.uorf_data: pd.DataFrame = None
		self.interorf_single: pd.DataFrame = None
		self.tmp_interorf_single_file: str = None
		self.interorfs_bed_df: pd.DataFrame = None
		self.tmp_interorfs_bed_sorted: str = None
		self.interorfs_dict: Dict[str, str] = None
		self.cds_dict: Dict[str, str] = None
		self.exons_gtf_df: pd.DataFrame = None

	def process_data(self, raw_data_processor: RawDataProcessor, uorf_processor: uORFDataProcessor) -> None:

		self.gtf_raw = raw_data_processor.gtf_raw
		self.fasta = raw_data_processor.fasta
		self.uorf_variation_data = uorf_processor.uorf_variation_data
		
		self.process_interorf_data()

	def process_interorf_data(self) -> None:
		self.exons_data, self.exons_data_file = self._get_exons_data(self.gtf_raw)
		self.cds_df = self._get_cds_gtf(self.gtf_raw)
		self.exons_gtf_dict = self._get_cds_list(self.cds_df)
		self.first_cds_df = self._get_first_cds(self.cds_df)
		self.uorf_data = self._get_uorf_data(self.uorf_variation_data, self.first_cds_df)
		self.interorf_single, self.tmp_interorf_single_file = self._get_interorf_data(self.uorf_data)
		self.interorfs_bed_df, self.tmp_interorfs_bed_sorted = self._intersect_interorf_with_exons(self.tmp_interorf_single_file, self.exons_data_file)
		self.interorfs_dict = self._get_interorfs_fasta(self.fasta, self.tmp_interorfs_bed_sorted)
		self.cds_dict = self._get_cds_regions(self.fasta, self.cds_df)
		self.exons_gtf_df = self._process_exons_gtf_dict(self.exons_gtf_dict, self.cds_dict)

	def _get_exons_data(self, gtf_raw: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
		exons_data = gtf_raw.loc[gtf_raw['type'] == 'exon']
		exons_data = exons_data.loc[:, ['chr', 'start', 'end', 'info', 'score', 'strand']]
		exons_data = exons_data.sort_values(by=['chr', 'start']).drop_duplicates()
		exons_data_file = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		exons_data.to_csv(exons_data_file.name, sep='\t', header=False, index=False)
		Logger.log_num_exons_data_records(exons_data.shape)
		return exons_data, exons_data_file
	
	def _get_cds_gtf(self, gtf_raw: pd.DataFrame) -> pd.DataFrame:
		cds_df = gtf_raw.loc[gtf_raw['type'] == 'CDS']
		cds_df = cds_df.sort_values(by=['chr', 'start']).drop_duplicates()
		Logger.log_num_cds_gtf_records(cds_df.shape)
		return cds_df
	
	def _get_cds_list(self, cds_df: pd.DataFrame) -> Dict[str, Dict[str, List[int]]]:
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
		Logger.log_num_cds_list_records(len(exons_gtf_dict))
		return exons_gtf_dict

	def _get_first_cds(self, cds_df: pd.DataFrame) -> pd.DataFrame:
		first_cds_df = cds_df.sort_values(['transcript', 'strand', 'start'])
		first_cds_df = first_cds_df.groupby('transcript').apply(lambda x: x.iloc[0] if x['strand'].iloc[0] == '+' else x.iloc[-1])
		first_cds_df.reset_index(drop=True, inplace=True)
		first_cds_df = first_cds_df.set_index('transcript')
		Logger.log_num_first_cds_records(first_cds_df.shape)
		return first_cds_df
	
	def _get_uorf_data(self, uorf_variation_data: pd.DataFrame, first_cds_df: pd.DataFrame) -> pd.DataFrame:
		uorf_data = uorf_variation_data.loc[:, ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_starts', 'name', 'name_and_trans']].copy()
		uorf_data.columns = ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', 'exons_sizes', 'exons_norm_starts', 'id', 'name_and_trans']
		uorf_data = uorf_data.map(lambda x: str(x) if isinstance(x, list) else x).drop_duplicates()
		uorf_data['exons_norm_starts'] = uorf_data['exons_norm_starts'].apply(eval)
		uorf_data['exons_sizes'] = uorf_data['exons_sizes'].apply(eval)
		uorf_data[['exons_starts', 'exons_norm_starts']] = uorf_data.apply(self._process_exons, axis=1, result_type='expand')
		uorf_data['exons_sizes'] = uorf_data.apply(lambda row: row['exons_sizes'] if row['strand'] == '+' else row['exons_sizes'][::-1], axis=1)
		uorf_data.index = uorf_data['id']
		uorf_data['cds_start'] = uorf_data['transcript'].apply(lambda x: self._get_true_cds_start(x, first_cds_df))
		Logger.log_num_uorf_data_records(uorf_data.shape)
		return uorf_data
	
	def _process_exons(self, row: pd.Series) -> Tuple[List[int], List[int]]:
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
		
	def _get_true_cds_start(self, tr_id: str, first_cds_df: pd.DataFrame) -> Optional[int]:
		if tr_id in first_cds_df.index:
			return first_cds_df.loc[tr_id, 'start'] if first_cds_df.loc[tr_id, 'strand'] == '+' else first_cds_df.loc[tr_id, 'end']
		return None
	
	def _get_interorf_data(self, uorf_data: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
		mask = ((uorf_data['strand'] == '+') & (uorf_data['orf_end'] < uorf_data['cds_start'])) | (uorf_data['cds_start'] < uorf_data['orf_start'])
		interorf_data = uorf_data[mask].copy()
		interorf_data.loc[interorf_data['strand'] == '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] == '+', ['orf_end', 'cds_start']].values
		interorf_data.loc[interorf_data['strand'] != '+', ['start', 'end']] = interorf_data.loc[interorf_data['strand'] != '+', ['cds_start', 'orf_start']].values
		interorf_single = interorf_data[['#CHROM', 'start', 'end', 'name_and_trans']].rename(columns={'#CHROM': 'chr'}).dropna()
		interorf_single[['start', 'end']] = interorf_single[['start', 'end']].astype(int)
		tmp_interorf_single_file = TemporaryFileManager.create('.bed', tmp_dir='/tmp')
		interorf_single.to_csv(tmp_interorf_single_file.name, sep='\t', header=False, index=False)
		Logger.log_num_interorf_single_records(interorf_single.shape)
		return interorf_single, tmp_interorf_single_file

	def _intersect_interorf_with_exons(self, tmp_interorf_single_file: str, exons_data_file: str) -> Tuple[pd.DataFrame, str]:

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
		Logger.log_num_interorfs_bed_records(interorfs_bed_df.shape)
		return interorfs_bed_df, tmp_interorfs_bed_sorted

	def _get_interorfs_fasta(self, fasta: str, tmp_interorfs_bed_sorted: str, strand: bool = True) -> Dict[str, str]:
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
		Logger.log_num_interorfs_fasta_records(len(interorfs_dict))
		return interorfs_dict

	def _get_cds_regions(self, fasta: str, cds_df: pd.DataFrame) -> Dict[str, str]:
		cds_bed_df = cds_df[['chr', 'start', 'end', 'strand']]
		cds_bed_df.loc[:, ['score', 'name']] = ['.', '.']
		cds_bed_df = cds_bed_df[['chr', 'start', 'end', 'name', 'score', 'strand']]
		cds_bed_df.loc[:, 'start'] = cds_bed_df['start'].astype(int, copy=False)
		cds_bed_df.loc[:, 'end'] = cds_bed_df['end'].astype(int, copy=False)
		cds_bed_df.loc[:, 'start'] -= 1
		cds_bed_df = cds_bed_df.sort_values(by=['chr', 'start'])
		cds_bed_df = cds_bed_df.drop_duplicates()

		cds_bed = TemporaryFileManager.create('.bed', tmp_dir='/tmp')

		cds_bed_df.to_csv(cds_bed.name, sep='\t', index=False, header=False)

		getfasta = Bedtools.getfasta(fasta, cds_bed.name, strand=True)

		cds_dict = SeqIO.to_dict(SeqIO.parse(getfasta, "fasta"))
		
		Logger.log_num_cds_regions_records(len(cds_dict))

		return cds_dict

	def _process_exons_gtf_dict(self, exons_gtf_dict: Dict[str, Dict[str, List[int]]], cds_dict: Dict[str, str]) -> pd.DataFrame:
		for k in exons_gtf_dict.keys():
			echrom = exons_gtf_dict[k]['chr']
			estrand = exons_gtf_dict[k]['strand']
			estarts = exons_gtf_dict[k]['exons_starts']
			esizes = exons_gtf_dict[k]['exons_sizes']
			exons_gtf_dict[k]['seq'] = str(self._create_complete_cds(echrom, estrand, estarts, esizes, cds_dict))
		exons_gtf_df = pd.DataFrame(exons_gtf_dict).T
		Logger.log_num_exons_gtf_records(exons_gtf_df.shape)
		return exons_gtf_df

	def _create_complete_cds(self, chrom: str, strand: str, eStarts: List[int], eSizes: List[int], cds_dict: Dict[str, str]) -> Optional[str]:
		if chrom in ['chrM', 'chrMT', 'M', 'MT']:
			return None
		out_seq = ''
		for i, j in zip(eStarts, eSizes):
			if strand == '+':
				out_seq += cds_dict[f'{chrom}:{i}-{i+j}({strand})'].seq
			else:
				out_seq += cds_dict[f'{chrom}:{i-j}-{i}({strand})'].seq
		return out_seq
