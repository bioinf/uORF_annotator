from typing import Dict, List, Optional, Tuple
from data_processor import RawDataProcessor, uORFDataProcessor, InterORFDataProcessor


from logger import Logger
from temporary_file_manager import TemporaryFileManager
from bedtools_wrapper import Bedtools
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import json
import re

import hashlib


class DataProcessor:
	def __init__(self, gtf_file_path, vcf_file_path, bed_file_path, fasta_file_path):
		self.raw_data_processor = RawDataProcessor()
		self.uorf_processor = uORFDataProcessor()
		self.interorf_processor = InterORFDataProcessor()

		self.gtf = gtf_file_path
		self.vcf = vcf_file_path
		self.bed = bed_file_path
		self.fasta = fasta_file_path

		self.bed_4col_info = '|utid|overlapping_type|dominance_type|codon_type'
		self.bed_4col_info_cols = self.bed_4col_info.split('|')[1:]

	def process_data(self):

		self._process_raw_data(self.vcf, self.bed, self.gtf, self.fasta)
		self._process_uorf_data(self.raw_data_processor, self.bed_4col_info_cols)
		self._process_interorf_data(self.raw_data_processor, self.uorf_processor)

	def _process_raw_data(self, vcf_file_path, bed_file_path, gtf_file_path, fasta_file_path):
		self.raw_data_processor.process_data(vcf_file_path, bed_file_path, gtf_file_path, fasta_file_path)

	def _process_uorf_data(self, raw_data_processor, bed_4col_info_cols):
		self.uorf_processor.process_data(bed_4col_info_cols, raw_data_processor)

	def _process_interorf_data(self, raw_data_processor, uorf_processor):
		self.interorf_processor.process_data(raw_data_processor, uorf_processor)