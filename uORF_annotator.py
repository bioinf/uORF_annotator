import subprocess as sp
from pathlib import Path
from tempfile import NamedTemporaryFile
import pandas as pd
import numpy as np
from math import ceil
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import click


def process(input_vcf, bed, fasta, bed_4col_info_cols, gtf, h) -> pd.core.frame.DataFrame:
	"""find main values for further annotation"""

	tmp_vcf1 = NamedTemporaryFile()
	sp.run(f'cut -f-8 {input_vcf} > {tmp_vcf1.name}', shell=True)

	# get intersection
	tmp_out = NamedTemporaryFile()
	sp.run(f'bedtools intersect -wo -a {tmp_vcf1.name} -b {bed} > {tmp_out.name}', shell=True)

	# preprocessing
	df = pd.read_table(tmp_out.name, header=None)
	df = df.dropna()

	df = df.loc[:, [0, 1, 3, 4, 7, 9, 10, 11, 13, 17, 18, 19]]
	df = df.rename({0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT',
					7: 'INFO', 9: 'orf_start', 10: 'orf_end',
					11: 'bed_anno', 13:'strand', 17:'n_exons', 
					18: 'exons_sizes', 19: 'exons_starts'}, axis=1)

	# string to list of int
	df['exons_sizes'] = df['exons_sizes'].apply(lambda x: [int(i) for i in x.rstrip(',').split(',')])
	df['exons_starts'] = df['exons_starts'].apply(lambda x: [int(i) for i in x.rstrip(',').split(',')])

	# get columns names for bed-4-field annotations
	bed_4col_values = df['bed_anno'].str.split('|', expand=True)
	bed_4col_values.columns = bed_4col_info_cols

	df = pd.concat([df, bed_4col_values], axis=1)

	# estimate distance
	df['dist_from_orf_to_snp'] = df.apply(lambda x: x.loc['POS'] - x.loc['orf_start'] - 1 if \
												x.loc['strand'] == '+' else \
												x.loc['orf_end'] - x.loc['POS'],
												axis = 1)

	# obtaining sequences from the corresponding reference genome
	tmp_fasta = NamedTemporaryFile()
	sp.run(f'rm -f {fasta}.fai && bedtools getfasta -s -fi {fasta} -bed {bed} > {tmp_fasta.name} 2>/dev/null',
			shell=True)
	record_dict = SeqIO.to_dict(SeqIO.parse(tmp_fasta.name, "fasta"))
	
	# df['INFO'] = df['INFO'].str.replace("[\'\"]$|^[\'\"]", "", regex=True)
	# annotate variants
	df[['symbol', 'consequence']] = df.apply(lambda x: annotate_variant(x, fasta_dict=record_dict), axis=1)
	
	# check overlapping with gtf annotation of CDS
	df['overlapped_CDS'] = ''
	if gtf is not None:
		df = check_overlapping(df, gtf, h)

	return df

def check_overlapping(df, gtf, h):
	"""mark whether variants intersect with the GTF-annotation"""

	# add rest obligate VCF-fields
	df = df.drop('overlapped_CDS', axis=1)
	check_cds_df = df.copy()
	check_cds_df['ID'] = '.'
	check_cds_df['QUAL'] = '.'
	check_cds_df['FILTER'] = '.'
	check_cds_df['FORMAT'] = '.'
	# check_cds_df['INFO'] = check_cds_df['INFO'].str.replace("[\'\"]$|^[\'\"]", "", regex=True)

	check_cds_df = check_cds_df.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

	# gtf to bed
	tmp_bed = NamedTemporaryFile()
	sp.run(f"""grep -v '#' {gtf} | awk '$3~/CDS/ {{OFS="\t";print$1,$4-1,$5,$3}}' > {tmp_bed.name}""", shell=True)

	# write VCF-header and VCF-body in temparary output file
	tmp_vcf2 = NamedTemporaryFile()
	tmp_tsv = NamedTemporaryFile()
	with open(tmp_vcf2.name, 'a') as w:
		w.write(h)
		check_cds_df.to_csv(w, sep='\t', index=None)

	# get intersection of current VCF and input GTF annotation
	sp.run(f'bedtools intersect -wo -a {tmp_vcf2.name} -b {tmp_bed.name} | cut -f1,2,4,5 | sort -k1,1 -k2,2n -u > {tmp_tsv.name}', shell=True)
	check_cds_df = pd.read_table(tmp_tsv.name, header=None)
	check_cds_df.columns = ['#CHROM', 'POS', 'REF', 'ALT']
	check_cds_df['overlapped_CDS'] = 'YES'
	df = pd.merge(df, check_cds_df, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left')
	df = df.fillna({'overlapped_CDS':'NO'})

	return df

def annotate_variant(x, fasta_dict) -> pd.core.series.Series:
	"""get a symbolic description of the mutations and consequences"""
	
	if x['ALT'] == '<DEL>':
		symb = '<DEL_L>'
		conseq = 'large_deletion'
		return pd.Series([symb, conseq], index=['symbol', 'consequence'])
	elif (len(x['REF']) > 1) & ((x['POS']<=x['orf_start']) | (x['POS']>x['orf_end'])):
		symb = '<DELB>'
		conseq = 'deletion_boundary'
		return pd.Series([symb, conseq], index=['symbol', 'consequence'])

	mnp_condition = (len(x['REF']) > 1) & (len(x['ALT']) > 1)
	if  mnp_condition:
		symb = '<MNP>'
		conseq = 'multinucleotide_polymorphism'
		return pd.Series([symb, conseq], index=['symbol', 'consequence'])

	# get fasta sequence of current ORF
	orf = fasta_dict[f"{x['#CHROM']}:{int(x.loc['orf_start'])}-{int(x.loc['orf_end'])}({x.loc['strand']})"]
	
	if x['n_exons'] == 1:
		# get position of nucl in codon
		x['codon_pos'] = x['dist_from_orf_to_snp'] % 3
		# get ORF-based coordinates
		begin = int(x['dist_from_orf_to_snp'] - x['codon_pos'])
		end = int(x['dist_from_orf_to_snp'] + (3 - x['codon_pos']))
		# get sequence of reference codon from ORF
		ref_codon = orf.seq[int(begin):int(end)]
		aa_pos = ceil(len(orf.seq)/3)

	# there is one or more introns
	else:   
		exons_ends = [i+j for i,j in zip(x['exons_sizes'],
										x['exons_starts'])]
		introns_starts = [int(i) for i in exons_ends][:-1] + [None]
		introns_ends = [int(i) for i in x['exons_starts'][1:]] + [None]

		block_ranges = pd.DataFrame({'ex_start': x['exons_starts'],
									'ex_end': exons_ends,
									'int_start': introns_starts,
									'int_end': introns_ends})
		if x['strand'] == '-':
			# swap start and end
			block_ranges_tr = pd.DataFrame({'ex_start': block_ranges['ex_end'],
											'ex_end': block_ranges['ex_start'],
											'int_start': block_ranges['int_end'],
											'int_end': block_ranges['int_start']})
			# flip dataframe
			block_ranges_tr = block_ranges_tr.iloc[::-1]
			block_ranges_tr = block_ranges_tr.reset_index()
			# subtract all coordinates from ORF length
			block_ranges_tr = block_ranges_tr.apply(lambda x: len(orf.seq) - x)
			# place NA from top to bottom as in initial data
			block_ranges_tr.loc[:, ['int_start', 'int_end']] = block_ranges_tr.loc[:, ['int_start', 'int_end']].shift(-1)
			# redifine data
			block_ranges = block_ranges_tr.copy()

		# add splice ranges
		block_ranges['splice_upstream_start'] = block_ranges['ex_start'] - 2
		block_ranges['splice_downstream_end'] = block_ranges['ex_end'] + 2

		# exclusive first and last exons
		block_ranges.loc[0, 'splice_upstream_start'] = None
		block_ranges.loc[len(block_ranges)-1, 'splice_downstream_end'] = None

		# check splice sites
		splice_condition = len(block_ranges[block_ranges.eval(f'((splice_upstream_start <= {x["dist_from_orf_to_snp"]}) & '
															f'({x["dist_from_orf_to_snp"]} < ex_start)) | '
															f'((ex_end <= {x["dist_from_orf_to_snp"]}) & '
															f'({x["dist_from_orf_to_snp"]} < splice_downstream_end))')]) != 0
		if splice_condition:
			symb = '<SP>'
			conseq = 'splice_variant'
			return pd.Series([symb, conseq], index=['symbol', 'consequence'])

		# is SNP in intron?
		intron_condition = len(block_ranges[block_ranges.eval(f'(int_start <= {x["dist_from_orf_to_snp"]}) & '
															f'({x["dist_from_orf_to_snp"]} < int_end)')]) != 0
		if intron_condition:
			symb = '<INT>'
			conseq = 'intronic_variant'
			return pd.Series([symb, conseq], index=['symbol', 'consequence'])

		else:
			# get intronic nucleotide counts before exon with SNP
			block_ranges['int_sum'] = (block_ranges['int_end'].shift(1) - block_ranges['int_start'].shift(1)).cumsum()                              
			intron_nucl_count = block_ranges[block_ranges.eval(f'(ex_start <= {x["dist_from_orf_to_snp"]}) & '
															f'({x["dist_from_orf_to_snp"]} < ex_end)')]['int_sum'].tolist()

			intron_nucl_count = 0 if (np.isnan(intron_nucl_count[0])) | (len(intron_nucl_count) == 0) else intron_nucl_count[0]

			# reformat snp distance from snp[orf_coord] to snp[cds_coord]
			x['dist_from_orf_to_snp'] = x['dist_from_orf_to_snp'] - intron_nucl_count

			# find codon position and coordinates of begin and end
			x['codon_pos'] = x['dist_from_orf_to_snp'] % 3
			begin = int(x['dist_from_orf_to_snp'] - x['codon_pos'])
			end = int(x['dist_from_orf_to_snp'] + (3 - x['codon_pos']))

			# get CDS sequence from ORF
			cds = Seq('')
			cds_list = [orf.seq[i['ex_start']:i['ex_end']] for _, i in block_ranges.loc[:, ['ex_start', 'ex_end']].iterrows()]
			for exon in cds_list:
				cds+=exon
			# get sequence of reference codon from cds
			ref_codon = cds[int(begin):int(end)]

			aa_pos = ceil(len(cds)/3)

	# get sequence of alternative codon
	ref_aa = str(ref_codon.translate())
	alt_codon = MutableSeq(ref_codon)

	# main condition
	diff = abs(len(x['REF']) - len(x['ALT']))
	# classify mutations into consequence groups
	if diff == 0:
		# insert mutation in reference codon and get new aminoacid or stop-codon
		if x['strand'] == '-':
			alt_codon[int(x['codon_pos'])] = str(Seq(x['ALT']).reverse_complement())
		else:
			alt_codon[int(x['codon_pos'])] = x['ALT']
		alt_codon = Seq(alt_codon)
		alt_aa = str(alt_codon.translate())

		# SNP variant
		if ref_aa == alt_aa:
			symb = f'{ref_aa}{aa_pos}{alt_aa}'
			conseq = 'synonymous_variant'
		elif (ref_aa != alt_aa) & (alt_aa != '*') & (ref_aa != '*'):
			symb = f'{ref_aa}{aa_pos}{alt_aa}'
			conseq = 'missense_variant'
		elif (ref_aa == '*') & (alt_aa != '*'):
			symb = f'*{aa_pos}{alt_aa}'
			conseq = 'stop_lost'
		elif (ref_aa != alt_aa) & (alt_aa == '*'):
			symb = f'{ref_aa}{aa_pos}*'
			conseq = 'stop_gained'
		elif (ref_codon in [Seq('ATG'), Seq('CTG')]) & (alt_codon not in [Seq('ATG'), Seq('CTG')]):
			symb = f'{ref_aa}{aa_pos}{alt_aa}'
			conseq = 'start_lost'

	# deletion/insertion/frameshift
	elif (diff % 3 == 0) & (len(x['REF']) > len(x['ALT'])):
		symb = '<DEL>'
		conseq = 'inframe_deletion'

	elif (diff % 3 == 0) & (len(x['REF']) < len(x['ALT'])):
		symb = '<INS>'
		conseq = 'inframe_insertion'

	elif diff % 3 != 0:
		symb = '<FS>'
		conseq = 'frameshift'

	# return two columns: symbolic designation of variants and consequence group
	return pd.Series([symb, conseq], index=['symbol', 'consequence'])


if __name__ == '__main__':

	@click.command()
	@click.option('--input_vcf', '-i', required=True)
	@click.option('--bed', '-b', required=True)
	@click.option('--fasta', '-f', required=True)
	@click.option('--gtf', '-g')
	@click.option('--output_tsv', '-ot')
	@click.option('--output_vcf', '-ov', required=True)
	def main(input_vcf, bed, fasta, gtf, output_tsv, output_vcf) -> None:

		# HARDCODED: annotations from 4 column of bed file
		bed_4col_info = '|utid|overlapping_type|codon_type'
		bed_4col_info_cols = bed_4col_info.split('|')[1:]

		# get VCF-header lines
		tmp_h = NamedTemporaryFile()
		sp.run(f"head -5000 {input_vcf} | grep -P '^##' > {tmp_h.name}", shell=True)
		with open(tmp_h.name) as f:
			h = f.read()
		# processing
		df = process(input_vcf, bed, fasta, bed_4col_info_cols, gtf, h)

		# write optional output tsv file
		if output_tsv is not None:
			df.to_csv(output_tsv, sep='\t', index=None)

		# get counts of variation types
		print(df['consequence'].value_counts())

		# write new CSQ like INFO line
		df['INFO_new'] = \
			[f'uBERT={int(x["orf_start"])}|'
			f'{int(x["orf_end"])}|'
			f'{x["strand"]}|'
			f'{x["symbol"]}|'
			f'{x["consequence"]}|'
			f'{x["overlapped_CDS"]}' for x in df.to_dict(orient='records')]

		# add fields from 4-field bed file
		df['INFO_new'] = df['INFO_new'].str.cat(df['bed_anno'], sep = '|')
		df['INFO_new'] = df['INFO_new'].str.replace(',uBERT=', ',')
		df['INFO_new'] = df['INFO_new'].astype(str) + ';'
		# set main VCF fields
		df = df.loc[:, ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'INFO_new']]

		# per variant multi-u-transcript annotation
		df = df.groupby(['#CHROM', 'POS', 'REF', 'ALT', 'INFO'], sort=False)['INFO_new'].apply(','.join)
		df = df.reset_index()

		# add rest obligate VCF-fields
		df['ID'] = '.'
		df['QUAL'] = '.'
		df['FILTER'] = '.'
		df['FORMAT'] = '.'
		df['INFO_new'] = df['INFO_new'].str.replace(',uBERT=', ',')
		# df['INFO'] = df['INFO'].str.replace("[\'\"]$|^[\'\"]", "", regex=True)
		df['INFO'] = df['INFO'].astype(str) + ';' + df['INFO_new'].astype(str)
		df = df.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]

		# remove output file before appending new data
		out_vcf = Path(output_vcf)
		if out_vcf.exists():
			out_vcf.unlink()

		# add header ##INFO uORF_annotator (uBERT) line
		h += \
		f'##INFO=<ID=uBERT,Number=.,Type=String,Description="Consequence uORF_annotator from uBERT. '
		f'Format: ORF_START|ORF_END|ORF_SYMB|ORF_CONSEQ|overlapped_CDS{bed_4col_info}">'
		# write VCF-header and VCF-body in output file
		with open(output_vcf, 'a') as w:
			w.write(h)
			df.to_csv(w, sep='\t', index=None)
	# run analysis
	main()
