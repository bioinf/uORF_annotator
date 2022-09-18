from email.policy import default
import subprocess as sp
from pathlib import Path
from tempfile import NamedTemporaryFile
import pandas as pd
import numpy as np
import re
import os
from collections import defaultdict
from math import ceil
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import click


def process(input_vcf, bed, fasta, bed_4col_info_cols, gtf, h, output_bed, atg_only) -> pd.core.frame.DataFrame:
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

        df['POS'] == df['POS'].astype(int)
        df['orf_start'] == df['orf_start'].astype(int)
        df['orf_end'] == df['orf_end'].astype(int)

        # estimate distance
        df['dist_from_orf_to_snp'] = df.apply(lambda x: x.loc['POS'] - x.loc['orf_start'] - 1 if \
                                                                                                x.loc['strand'] == '+' else \
                                                                                                x.loc['orf_end'] - x.loc['POS'],
                                                                                                axis = 1)
        df['dist_from_orf_to_snp'] = df['dist_from_orf_to_snp'].astype(int)

        # obtaining sequences from the corresponding reference genome
        tmp_fasta = NamedTemporaryFile()
        sp.run(f'rm -f {fasta}.fai && bedtools getfasta -s -fi {fasta} -bed {bed} > {tmp_fasta.name} 2>/dev/null',
                        shell=True)
        uorf_dict = SeqIO.to_dict(SeqIO.parse(tmp_fasta.name, "fasta"))

        tmp = NamedTemporaryFile()
        tmp1 = NamedTemporaryFile()
        tmp_exons_bed = NamedTemporaryFile()
        tmp_exons_gtf = NamedTemporaryFile()
        tmp_exons_sorted = NamedTemporaryFile()
        tmp_interorf_single = NamedTemporaryFile()
        tmp_io_exon_isec_tab = NamedTemporaryFile()
        tmp_interorfs_bed = NamedTemporaryFile()
        tmp_interorfs_bed_sorted = NamedTemporaryFile()
        tmp_interorfs_full_seq = NamedTemporaryFile()
        tmp_interorfs_fasta = NamedTemporaryFile()
        closest_cds_data = NamedTemporaryFile()

        df['zero'] = 0
        df = df.astype({'orf_start':'int32', 'orf_end':'int32'})
        df['name'] = df['#CHROM'] + ':' + df['orf_start'].astype(str) + '-' + df['orf_end'].astype(str) + '(' + df['strand'] + ')'
        df['transcript'] = [re.findall('[NMENST]+_*\d+', x)[0] for x in df['bed_anno']]
        df['name_and_trans'] = [f'{x}/{y}' for x, y in zip(df['name'], df['transcript'])]
        df.to_csv('check_temp_df.tsv', header=None, sep='\t', index=False)
        df.loc[:, ['#CHROM', 'orf_start', 'orf_end', 'name_and_trans', 'zero', 'strand']]\
                .to_csv(tmp.name, header=None, sep='\t', index=False)

        sp.run(f'sort -k1,1 -k2,2n {tmp.name} > {tmp1.name}', shell=True)
        
        source_gtf = pd.read_csv(gtf, sep='\t', comment='#', header=None)
        source_gtf.columns = ['chr', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'info']
        source_gtf['start'] = [int(x) - 1 for x in source_gtf['start']]
        source_gtf['end'] = [int(x) for x in source_gtf['end']]
        ok_rows = [bool(re.search('(NM_|ENST)\d+', x)) for x in source_gtf['info']]
        source_gtf = source_gtf.loc[ok_rows]
        source_gtf['transcript'] = [re.findall('transcript_id \"((NM_|ENST)\d+)', x)[0][0] for x in source_gtf['info']]

        exons_gtf = source_gtf.loc[source_gtf['type'] == 'exon']
        export_exons = exons_gtf.loc[:, ['chr', 'start', 'end', 'transcript', 'score', 'strand']]
        export_exons.to_csv(tmp_exons_bed.name, sep='\t', header=False, index=False)
        sp.run(f'sort -k1,1 -k4,4n {tmp_exons_bed.name} | uniq > {tmp_exons_sorted.name}', shell=True)

        cds_gtf = source_gtf.loc[source_gtf['type'] == 'CDS']

        uorf_df = df.loc[:, ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', \
                        'exons_sizes', 'exons_starts', 'name', 'name_and_trans']]
        uorf_df.columns = ['#CHROM', 'transcript', 'strand', 'orf_start', 'orf_end', \
                        'exons_sizes', 'exons_norm_starts', 'id', 'name_and_trans']
        uorf_df = uorf_df.loc[uorf_df.astype(str).drop_duplicates().index]
        uorf_exons = []
        uorf_norm_exons = []
        for idx, row_data in uorf_df.iterrows():
                out_starts = []
                out_norm_starts = []
                orf_len = row_data.loc['orf_end'] - row_data.loc['orf_start']
                if row_data.loc['strand'] == '+':
                        for eStart in row_data.loc['exons_norm_starts']:
                                out_starts.append(row_data.loc['orf_start'] + eStart)
                                out_norm_starts.append(eStart)
                else:
                        for eStart, eSize in zip(row_data.loc['exons_norm_starts'], row_data['exons_sizes']):
                                out_starts.append(row_data.loc['orf_start'] + eStart + eSize)
                                out_norm_starts.append(orf_len - eStart - eSize)
                        out_starts = list(reversed(out_starts))
                        out_norm_starts = list(reversed(out_norm_starts))
                uorf_exons.append(list(out_starts))
                uorf_norm_exons.append(list(out_norm_starts))
        uorf_df['exons_starts'] = uorf_exons
        uorf_df['exons_norm_starts'] = uorf_norm_exons
        uorf_esizes = []
        for idx, row_data in uorf_df.iterrows():
                if row_data['strand'] == '+':
                        uorf_esizes.append(list(row_data['exons_sizes']))
                else:
                        uorf_esizes.append(list(reversed(row_data['exons_sizes'])))
        uorf_df['exons_sizes'] = uorf_esizes

        #uorf_df['exons_starts'] = [(o + x for x in e) for o, e in zip(uorf_df['orf_start'], uorf_df['exons_norm_starts'])]
        uorf_df.index = uorf_df['id']
        #print(uorf_df.loc['chr11:31806911-31810840(-)'])

        first_cds = {}
        for idx, row_data in cds_gtf.iterrows():
                if row_data['transcript'] not in first_cds:
                        first_cds[row_data['transcript']] = row_data
                        continue
                if (row_data['strand'] == '-') and (first_cds[row_data['transcript']]['start'] < row_data['start']):
                     first_cds[row_data['transcript']] = row_data
        first_cds_df = pd.DataFrame(first_cds).T
        #print(first_cds_df)
        first_cds_df.index = first_cds_df['transcript']

        def get_true_cds_start(tr_id, first_cds_df=first_cds_df):
                if tr_id not in first_cds_df.index:
                        return None
                gtf_row = first_cds_df.loc[tr_id]
                if gtf_row['strand'] == '+':
                        return gtf_row['start']
                return gtf_row['end']

        uorf_df['cds_start'] = [get_true_cds_start(x) for x in uorf_df['transcript']]
        #print(uorf_df)

        interorf_data = []
        for idx, row_data in uorf_df.iterrows():
                if (row_data['strand'] == '+') and (row_data['orf_end'] < row_data['cds_start']):
                        interorf_data.append(row_data[['#CHROM', 'orf_end', 'cds_start', 'name_and_trans']].values)
                elif row_data['cds_start'] < row_data['orf_start']:
                        interorf_data.append(row_data[['#CHROM', 'cds_start', 'orf_start', 'name_and_trans']].values)
        interorf_single = pd.DataFrame(interorf_data)
        interorf_single = interorf_single.dropna()
        interorf_single.columns = ['chr', 'start', 'end', 'name_and_trans']
        interorf_single['start'] = [int(x) for x in interorf_single['start']]
        interorf_single['end'] = [int(x) for x in interorf_single['end']]
        interorf_single.to_csv(tmp_interorf_single.name, sep='\t', header=False, index=False)
        
        cmd = f"""bedtools intersect -wo -a {tmp_interorf_single.name} -b {tmp_exons_sorted.name} | tee test_isec.bed > {tmp_io_exon_isec_tab.name}"""
        sp.run(cmd, shell=True)

        with open(tmp_io_exon_isec_tab.name, 'r') as isec_handle, open(tmp_interorfs_bed.name, 'w') as isec_out:
                for line in isec_handle:
                        content = line.strip().split('\t')
                        uorf_tr_id = content[3].split('/')[1]
                        if uorf_tr_id == content[7]:
                                rstart = max(int(content[1]), int(content[5]))
                                rend = min(int(content[2]), int(content[6]))
                                out_line = f'{content[0]}\t{rstart}\t{rend}\t{content[3]}\t{content[8]}\t{content[9]}\n'
                                isec_out.write(out_line)
        
        command = f"""cat {tmp_interorfs_bed.name} | sort -k1,1 -k2,2n | uniq | tee test_io_exons.bed > {tmp_interorfs_bed_sorted.name}"""
        sp.run(command, shell=True)

        command = f"""rm -f {fasta}.fai && bedtools getfasta -bed {tmp_interorfs_bed_sorted.name} -fi {fasta} -s -name > {tmp_interorfs_fasta.name}""" #| sed 's/::.*//g' > {tmp_interorfs_fasta.name}"""
        sp.run(command, shell=True)
        

        ### INTERORF DATA ###
        interorfs_bed_list = [line.rstrip().split('\t') for line in open(tmp_interorfs_bed_sorted.name).read().split('\n') if line != '' ]
        interorfs_bed_dict = defaultdict(lambda: defaultdict(list))
        for bed_line in interorfs_bed_list:
                key = bed_line[3]
                interorfs_bed_dict[key]['chr'] = bed_line[0]
                interorfs_bed_dict[key]['strand'] = bed_line[-1]
                interorfs_bed_dict[key]['exons_sizes'].append(int(bed_line[2])-int(bed_line[1]))
                if interorfs_bed_dict[key]['strand'] ==  '+':
                        interorfs_bed_dict[key]['exons_starts'].append(int(bed_line[1]))
                elif interorfs_bed_dict[key]['strand'] ==  '-':
                        interorfs_bed_dict[key]['exons_starts'].append(int(bed_line[2]))
        for k in interorfs_bed_dict.keys():
                if interorfs_bed_dict[k]['strand'] == '+':
                        starts_sorted = sorted(interorfs_bed_dict[k]['exons_starts'])
                        interorfs_bed_dict[k]['exons_sizes'] = [x for _, x in sorted(zip(starts_sorted, interorfs_bed_dict[k]['exons_sizes']))]
                        contig = interorfs_bed_dict[k]['chr']
                        interorfs_bed_dict[k]['exons_norm_starts'] = [i-starts_sorted[0] for i in starts_sorted]
                elif interorfs_bed_dict[k]['strand'] == '-':
                        size_dict = {k: v for k, v in zip(interorfs_bed_dict[k]['exons_starts'], interorfs_bed_dict[k]['exons_sizes'])}
                        starts_sorted = sorted(interorfs_bed_dict[k]['exons_starts'], reverse=True)
                        interorfs_bed_dict[k]['exons_sizes'] = [size_dict[x] for x in starts_sorted]
                        contig = interorfs_bed_dict[k]['chr']
                        interorfs_bed_dict[k]['exons_norm_starts'] = [abs(i-starts_sorted[0]) for i in starts_sorted]
                        interorfs_bed_dict[k]['exons_starts'] = starts_sorted
        interorfs_bed_df = pd.DataFrame(interorfs_bed_dict).T
        # interorfs_bed_df.columns = ['contig', 'start', 'end']
        interorfs_bed_df['id'] = interorfs_bed_df.index
        #print(' ======= interORF DATA ======= ')
        #print(interorfs_bed_df)
        #print(interorfs_bed_df.loc['chr11:31810829-31811221(-)/NM_001604'])
        #####################

        seen_parts = set()
        records = []

        ############ BIG CHUNK OF NEW CODE ####################################

        interorfs_full_seq = defaultdict(str)

        with open(tmp_interorfs_fasta.name, "r") as split_fasta, open(tmp_interorfs_full_seq.name, 'w') as full_fasta:
            for line in split_fasta:
                if line.startswith('>'):
                    record_tmp_name = line[1:]
                    uorf_tmp_name = record_tmp_name.split('::')[0]
                    if record_tmp_name in seen_parts:
                        skip = True
                    else:
                        seen_parts.add(record_tmp_name)
                        skip = False
                else:
                    if skip:
                        continue
                    if '(-)' in uorf_tmp_name:
                        interorfs_full_seq[uorf_tmp_name] = line.strip() + interorfs_full_seq[uorf_tmp_name]
                    else:
                        interorfs_full_seq[uorf_tmp_name] = interorfs_full_seq[uorf_tmp_name] + line.strip()
            for interorf_seq in interorfs_full_seq:
                full_fasta.write(f">{interorf_seq}\n{interorfs_full_seq[interorf_seq]}\n")

        ############# END OF CHUNK ############################################


        cmd = f'grep -P "\tCDS\t" {gtf} | sort -k1,1 -k4,4n | uniq - > {tmp_exons_gtf.name}'
        sp.run(cmd, shell=True)

        # make per transcript exon dict (CDS-data)
        exons_gtf_list = [line.rstrip().split('\t') for line in open(tmp_exons_gtf.name).read().split('\n') if line.startswith('chr')]
        exons_gtf_dict = defaultdict(lambda: defaultdict(list))
        idx = [False for _ in exons_gtf_list[0]]
        idx[0], idx[3], idx[4], idx[6] = True, True, True, True
        for gtf_line in exons_gtf_list:
                key = re.search('transcript_id \"([^\.]+)', gtf_line[8]).group(1).strip('"')
                values = [gtf_line_sub for i, gtf_line_sub in zip(idx, gtf_line) if i]
                exons_gtf_dict[key]['chr'] = values[0]
                exons_gtf_dict[key]['strand'] = values[-1]
                exons_gtf_dict[key]['exons_sizes'].append(int(values[2])-int(values[1]) + 1)
                if exons_gtf_dict[key]['strand'] == '+':
                        exons_gtf_dict[key]['exons_starts'].append(int(values[1]) - 1)
                elif exons_gtf_dict[key]['strand'] == '-':
                        exons_gtf_dict[key]['exons_starts'].append(int(values[2]))
        for k in exons_gtf_dict.keys():
                if exons_gtf_dict[k]['strand'] == '+':
                        starts_sorted = sorted(exons_gtf_dict[k]['exons_starts'])
                        exons_gtf_dict[k]['exons_sizes'] = [x for _, x in sorted(zip(starts_sorted, exons_gtf_dict[k]['exons_sizes']))]
                        contig = exons_gtf_dict[k]['chr']
                        exons_gtf_dict[k]['id'] = f"{contig}:{starts_sorted[0]}(+)"
                        exons_gtf_dict[k]['exons_norm_starts'] = [i-starts_sorted[0] for i in starts_sorted]
                elif exons_gtf_dict[k]['strand'] == '-':
                        size_dict = {k: v for k, v in zip(exons_gtf_dict[k]['exons_starts'], exons_gtf_dict[k]['exons_sizes'])}
                        starts_sorted = sorted(exons_gtf_dict[k]['exons_starts'], reverse=True)
                        exons_gtf_dict[k]['exons_starts'] = starts_sorted
                        exons_gtf_dict[k]['exons_sizes'] = [size_dict[x] for x in starts_sorted]
                        contig = exons_gtf_dict[k]['chr']
                        exons_gtf_dict[k]['id'] = f"{contig}:{starts_sorted[0]}(-)"
                        exons_gtf_dict[k]['exons_norm_starts'] = [abs(i-starts_sorted[0]) for i in starts_sorted]

        ### Get sequences of CDS
        cds_fasta = NamedTemporaryFile()
        get_cds_command = f"""rm -f {fasta}.fai && cat {gtf} | awk -F '\t' '$3~/CDS/' | awk -F '\t' '{{OFS="\t";print$1,$4-1,$5,$3,0,$7}}' | sort -k1,1 -k2,2n | uniq | bedtools getfasta -bed - -fi {fasta} -s > {cds_fasta.name}"""
        sp.run(get_cds_command, shell=True)
        cds_seqs = []
        for record in SeqIO.parse(cds_fasta.name, "fasta"):  
            cds_seqs.append(record)

        cds_dict = SeqIO.to_dict(cds_seqs)
        
        def create_complete_cds(chrom, strand, eStarts, eSizes):
                if chrom in ['chrM', 'chrMT', 'M', 'MT']:
                        return None
                out_seq = ''
                for i, j in zip(eStarts, eSizes):
                        if strand == '+':
                                out_seq += cds_dict[f'{chrom}:{i}-{i+j}({strand})'].seq
                        else:
                                out_seq += cds_dict[f'{chrom}:{i-j}-{i}({strand})'].seq
                return out_seq
        
        for k in exons_gtf_dict.keys():
                echrom = exons_gtf_dict[k]['chr']
                estrand = exons_gtf_dict[k]['strand']
                estarts = exons_gtf_dict[k]['exons_starts']
                esizes = exons_gtf_dict[k]['exons_sizes']
                exons_gtf_dict[k]['seq'] = str(create_complete_cds(echrom, estrand, estarts, esizes))
        
        closest_cds_dict = defaultdict(str)
        with open(closest_cds_data.name, 'r') as uorf_cds_handle:
                for line in uorf_cds_handle:
                        content = line.strip().split('\t')
                        if content[5] == "+":
                                exStart = content[7]
                        else:
                                exStart = content[8]
                        closest_cds_dict[content[3]] = f"{content[0]}:{exStart}({content[5]})"

        exons_gtf_df = pd.DataFrame(exons_gtf_dict).T
        #exons_gtf_df.index = exons_gtf_df['id']
        #print(' ======= CDS DATA ======= ')
        #print(exons_gtf_df)
        exons_gtf_df.to_csv('look_at_me.tsv', sep='\t', index=False)
        ################

                

        for record in SeqIO.parse(tmp_interorfs_full_seq.name, "fasta"):  
            records.append(record)

        interorfs_dict = SeqIO.to_dict(records)

        # annotate variants
        df[['symbol', 'consequence', 'main_cds_effect']] = df.apply(lambda x: annotate_variant(x, \
                                                                                        uorf_dict=uorf_dict, \
                                                                                        interorfs_dict=interorfs_dict, \
                                                                                        interorfs_bed_df=interorfs_bed_df, \
                                                                                        uorf_df=uorf_df, \
                                                                                        cds_df=exons_gtf_df, \
                                                                                        closest_cds_dict=closest_cds_dict, \
                                                                                        output_bed=output_bed, atg_only=atg_only),
                                                                                        axis=1)

        # check overlapping with gtf annotation of CDS
        df['overlapped_CDS'] = ''
        df['additional_info'] = ''
        if gtf is not None:
                df = check_overlapping(df, gtf, h, fasta, bed_4col_info_cols)

        return df

def check_overlapping(df, gtf, h, fasta, bed_4col_info_cols) -> pd.core.frame.DataFrame:
        """mark whether variants intersect with the GTF-annotation"""

        # add rest obligate VCF-fields
        df = df.drop('overlapped_CDS', axis=1)
        check_cds_df = df.copy()
        check_cds_df['ID'] = '.'
        check_cds_df['QUAL'] = '.'
        check_cds_df['FILTER'] = '.'
        check_cds_df['FORMAT'] = '.'

        check_cds_df = check_cds_df.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

        # convert gtf to bed
        # read gtf as tsv, skip all lines starting with '#'
        df_gtf = pd.read_table(gtf, comment='#', header=None)
        # subset only CDS annotation
        df_gtf = df_gtf[df_gtf[2] == 'CDS']
        # extract gene_id without quotes
        df_gtf[8] = df_gtf[8].str.extract('gene_id [\'|\"]?([^;\"\']+)')

        # create dict {'<gene_id>': [<chr>, <start>, <end>, info>]}
        d = defaultdict(list)
        for _, i in df_gtf.iterrows():
                if i[8] in d:
                        if d[i[8]][1] > i[3]-1:
                                d[i[8]][1] = i[3]-1
                        if d[i[8]][2] < i[4]:
                                d[i[8]][2] = i[4]
                        continue
                else:
                        d[i[8]] = [i[0], i[3]-1, i[4], i[6]]
        # make dataframe from dict
        d = pd.DataFrame(d).T
        d[3] = d.index
        # write temporary bed file
        tmp_bed1 = NamedTemporaryFile()
        pd.DataFrame(d).to_csv(tmp_bed1.name, sep='\t', index=None, header=None)
        # write VCF-header and VCF-body in temparary output file
        tmp_vcf2 = NamedTemporaryFile()
        tmp_tsv = NamedTemporaryFile()
        with open(tmp_vcf2.name, 'a') as w:
                w.write(h)
                check_cds_df.to_csv(w, sep='\t', index=None)

        # get intersection of current VCF and input GTF annotation and obtain information about ORF interactions
        tmp_bed2 = NamedTemporaryFile()
        sp.run(f'sort -k1,1 -k2,2n {tmp_bed1.name} | bedtools merge -i - | sort -k1,1 -k2,2n > {tmp_bed2.name}', shell=True)
        sp.run(f'bedtools intersect -wo -a {tmp_vcf2.name} -b {tmp_bed2.name} | cut -f1,2,4,5 | sort -k1,1 -k2,2n > {tmp_tsv.name}', shell=True)
        check_cds_df = pd.read_table(tmp_tsv.name, header=None)
        check_cds_df.columns = ['#CHROM', 'POS', 'REF', 'ALT']
        df.to_csv('df.tsv', sep='\t', index=None)
        check_cds_df['in_known_ORF'] = 'YES'
        check_cds_df.to_csv('check_cds_df.tsv', sep='\t', index=None)

        df = pd.merge(df, check_cds_df, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left')
        df = df.fillna({'in_known_ORF':'NO'})
        del check_cds_df

        # get intersection of current VCF and input GTF annotation and obtain information about CDS interactions
        tmp_bed3 = NamedTemporaryFile()
        tmp_tsv2 = NamedTemporaryFile()
        sp.run(f"""grep -v '#' {gtf} | awk '$3~/CDS/ {{OFS="\t";print$1,$4-3,$5+2,$3}}' | sort -k1,1 -k2,2n | bedtools merge -i - > {tmp_bed3.name}""", shell=True)
        sp.run(f'bedtools intersect -wo -a {tmp_vcf2.name} -b {tmp_bed3.name} | cut -f1,2,4,5 | sort -k1,1 -k2,2n | uniq > {tmp_tsv2.name}', shell=True)
        check_cds_df2 = pd.read_table(tmp_tsv2.name, header=None)
        check_cds_df2.columns = ['#CHROM', 'POS', 'REF', 'ALT']
        df.to_csv('df.tsv', sep='\t', index=None)
        check_cds_df2['in_known_CDS'] = 'YES'
        check_cds_df2.to_csv('check_cds_df.tsv', sep='\t', index=None)
        df = pd.merge(df, check_cds_df2, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left')
        df = df.fillna({'in_known_CDS':'NO'})
        df.loc[df['in_known_ORF'] == 'NO', 'in_known_CDS'] = 'NO'

        return df


def annotate_variant(x, uorf_dict, interorfs_dict, interorfs_bed_df, uorf_df, cds_df, closest_cds_dict, output_bed, atg_only) -> pd.core.series.Series:
        """get a symbolic description of the mutations and consequences"""

        main_cds_effect = ''

        # checking for the effect of large deletions
        if x['ALT'] == '<DEL>':
                symb = '<DEL_L>'
                conseq = 'large_deletion'
                return pd.Series([symb, conseq, main_cds_effect], index=['symbol', 'consequence', 'main_cds_effect'])

        elif (len(x['REF']) > 1) & ((x['POS']<=x['orf_start']) | (x['POS']>x['orf_end'])):
                symb = '<DELB>'
                conseq = 'deletion_boundary'
                return pd.Series([symb, conseq, main_cds_effect], index=['symbol', 'consequence', 'main_cds_effect'])

        # excluding of multinucleotide polymorphism
        mnp_condition = (len(x['REF']) > 1) & (len(x['ALT']) > 1)
        if  mnp_condition:
                symb = '<MNP>'
                conseq = 'multinucleotide_polymorphism'
                return pd.Series([symb, conseq, main_cds_effect], index=['symbol', 'consequence', 'main_cds_effect'])

        # get fasta sequence of current ORF
        orf = uorf_dict[f"{x['#CHROM']}:{int(x.loc['orf_start'])}-{int(x.loc['orf_end'])}({x.loc['strand']})"]
        cds = Seq(orf.seq)

        diff = abs(len(x['REF']) - len(x['ALT']))
        signed_diff = len(x['REF']) - len(x['ALT'])
        del_end_intron = False

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
                        return pd.Series([symb, conseq, main_cds_effect], index=['symbol', 'consequence', 'main_cds_effect'])

                # is SNP in intron?
                intron_condition = len(block_ranges[block_ranges.eval(f'(int_start <= {x["dist_from_orf_to_snp"]}) & '
                                                                                                                        f'({x["dist_from_orf_to_snp"]} < int_end)')]) != 0
                if intron_condition:
                        symb = '<INT>'
                        conseq = 'intronic_variant'
                        return pd.Series([symb, conseq, main_cds_effect], index=['symbol', 'consequence', 'main_cds_effect'])

                else:
                        # get intronic nucleotide counts before exon with SNP
                        block_ranges['int_sum'] = (block_ranges['int_end'].shift(1) - block_ranges['int_start'].shift(1)).cumsum()                              
                        intron_nucl_count = block_ranges[block_ranges.eval(f'(ex_start <= {x["dist_from_orf_to_snp"]}) & '
                                                                                                                        f'({x["dist_from_orf_to_snp"]} < ex_end)')]['int_sum'].tolist()

                        intron_nucl_count = 0 if (np.isnan(intron_nucl_count[0])) | (len(intron_nucl_count) == 0) else intron_nucl_count[0]

                        if len(x['REF']) > len(x['ALT']):
                                if x['strand'] == "+":
                                        del_end_pos = x["dist_from_orf_to_snp"] + diff
                                else:
                                        del_end_pos = x["dist_from_orf_to_snp"] - diff
                                        del_end_intron = len(block_ranges[block_ranges.eval(f'(int_start <= {del_end_pos}) & '
                                                                                                                        f'({del_end_pos} < int_end)')]) != 0
                

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
                        alt_seq = cds[:-3] + str(alt_codon)
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
                if del_end_intron:
                        symb = '<EMD>'
                        conseq = 'deletion_at_exon_margin'
                else:   
                        symb = '<FS>'
                        conseq = 'frameshift'
                        x['dist_from_orf_to_snp'] = int(x['dist_from_orf_to_snp'])
                        if len(x['ALT']) > len(x['REF']):
                                if x['strand'] == "+":
                                        alt_seq = str(cds[:x['dist_from_orf_to_snp']] + Seq(x['ALT']) + cds[x['dist_from_orf_to_snp'] + 1:])
                                elif x['strand'] == '-':
                                        alt_seq = str(cds[:x['dist_from_orf_to_snp']] + Seq(x['ALT']).reverse_complement() + cds[x['dist_from_orf_to_snp'] + 1:])

                                alt_seq = Seq(alt_seq)
                        else:
                                if x['strand'] == "+":
                                        alt_seq = str(cds[:x['dist_from_orf_to_snp']] + Seq(x['ALT']) + cds[x['dist_from_orf_to_snp'] + len(x['ALT']):])
                                elif x['strand'] == '-':
                                        alt_seq = str(cds[:x['dist_from_orf_to_snp'] - diff] + Seq(x['ALT']).reverse_complement() + cds[x['dist_from_orf_to_snp'] + 1:])

        if (conseq == 'stop_lost') | (conseq == 'frameshift'):

                def get_intron_ln(exon_starts, exon_sizes):
                        if len(exon_sizes) == 1:
                                return []
                        return [exon_starts[i + 1] - exon_starts[i] - exon_sizes[i] for i in range(len(exon_sizes) - 1)]

                
                uorf_name = f"{x['#CHROM']}:{int(x.loc['orf_start'])}-{int(x.loc['orf_end'])}({x.loc['strand']})"
                uorf_len = int(x.loc['orf_end']) - int(x.loc['orf_start'])

                uorf_df_row = uorf_df.drop_duplicates(subset='id', keep='last').loc[uorf_name].copy()
                ORFExonSizes = uorf_df_row['exons_sizes'].copy()
                ORFExonStarts = uorf_df_row['exons_starts'].copy()
                ORFExonNormStarts = uorf_df_row['exons_norm_starts'].copy()

                if x['strand'] == '-':
                        #ORFExonNormStarts = list(reversed([uorf_len - x - y for x, y in zip(ORFExonSizes, ORFExonNormStarts)]))
                        #ORFExonSizes = list(reversed(ORFExonSizes))
                        true_uorf_start = x.loc['orf_end']
                        true_uorf_end = x.loc['orf_start']
                else:
                        true_uorf_start = x.loc['orf_start']
                        true_uorf_end = x.loc['orf_end']

                ORFIntronSizes = get_intron_ln(exon_starts=ORFExonNormStarts, exon_sizes=ORFExonSizes)
                
                transcript_id = uorf_df.loc[uorf_name]['transcript']
                try:
                        closest_cds_row = cds_df.loc[transcript_id].copy()
                        CDSExonSizes = closest_cds_row.loc['exons_sizes'].copy()
                        CDSExonStarts = closest_cds_row.loc['exons_starts'].copy()
                        CDSExonNormStarts = closest_cds_row.loc['exons_norm_starts'].copy()
                        CDSIntronSizes = get_intron_ln(exon_starts=CDSExonNormStarts, exon_sizes=CDSExonSizes)
                        CDSSeq = closest_cds_row.loc['seq']
                except KeyError:
                        if x.loc['overlapping_type'] == 'non-overlapping':
                                print(f'WARN Closest CDS is not the first one for non-overlappiong uORF {uorf_name}!')
                        return pd.Series([symb, conseq, 'unassigned'], index=['symbol', 'consequence', 'main_cds_effect'])

                try:
                        io_name = x.loc["name_and_trans"]
                        _, _, IOExonSizes, IOExonStarts, IOExonNormStarts, _ = interorfs_bed_df.loc[io_name].copy()
                        interorf = interorfs_dict[io_name].seq
                except KeyError:
                        IOExonSizes, IOExonStarts, IOExonNormStarts = [[], [], []]
                        interorf = ''
                
                if x.loc['strand'] == '+':
                        IOLen = CDSExonStarts[0] - true_uorf_end
                else:
                        IOLen = true_uorf_end - CDSExonStarts[0]
                if IOLen <= 0:
                        if x.loc['overlapping_type'] == 'non-overlapping':
                                print('WARN Saw the end of non-overlapping uORF inside CDS')
                        remaining_exons = CDSExonSizes.copy()
                        remaining_introns = CDSIntronSizes.copy()
                        for eID, eStart in enumerate(CDSExonStarts):
                                eSize = CDSExonSizes[eID]
                                if x.loc['strand'] == '+':
                                        eEnd = eStart + eSize
                                        pos_cond = (eEnd < true_uorf_end)
                                else:
                                        eEnd = eStart - eSize
                                        pos_cond = (eEnd >= true_uorf_end)
                                if pos_cond:
                                        remaining_exons = remaining_exons[1:]
                                        remaining_introns = remaining_introns[1:]
                                        CDSSeq = CDSSeq[eSize:]
                                else:
                                        extension = abs(eEnd - true_uorf_end)
                                        remove_piece = (remaining_exons[0] - extension)
                                        ORFExonSizes[-1] += extension 
                                        CDSExonSizes = remaining_exons[1:]
                                        CDSIntronSizes = remaining_introns
                                        CDSSeq = CDSSeq[remove_piece:]
                                        #if len(remaining_introns) > 1:
                                        #        CDSIntronSizes = remaining_introns[1:]
                                        #else:
                                        #        CDSIntronSizes = []
                                        break
                        AllExons = ORFExonSizes + CDSExonSizes
                        AllIntrons = ORFIntronSizes + CDSIntronSizes

                if (len(IOExonSizes) == 0) and (IOLen > 0):
                        AllExons = ORFExonSizes + CDSExonSizes
                        AllIntrons = ORFIntronSizes + [IOLen] + CDSIntronSizes
                elif IOLen > 0:
                        # Getting true end of the inter-ORF region
                        if x.loc['strand'] == '+':
                                true_io_end = IOExonStarts[-1] + IOExonSizes[-1]
                        else:
                                true_io_end = IOExonStarts[-1] - IOExonSizes[-1]
                        
                        # Getting IO-specific introns - left-most and right-most
                        if len(IOExonSizes) > 1:
                                IOInternalIntrons = get_intron_ln(exon_starts=IOExonNormStarts, exon_sizes=IOExonSizes)
                        else:
                                IOInternalIntrons = []                        
                        left_intron_size = abs(IOExonStarts[0] - true_uorf_end)
                        right_intron_size = abs(true_io_end - CDSExonStarts[0])

                        #if io_name == 'chr11:31806911-31810840(-)/NM_000280':
                        #        print(uorf_name)
                        #        print(uorf_df.loc[uorf_name])
                        #        print(list(ORFExonStarts))
                        #        print(ORFExonNormStarts)
                        #        print(ORFExonSizes)
                        #        print(IOExonStarts)
                        #        print(IOExonNormStarts)
                        #        print(IOExonSizes)
                        #        print(CDSExonStarts)
                        #        print(CDSExonNormStarts)
                        #        print(IOInternalIntrons)
                        #        print('## #INTRONS')
                        #        print(ORFIntronSizes)
                        #        print(left_intron_size)
                        #        print(IOInternalIntrons)
                        #        print(right_intron_size)
                        #        print(CDSIntronSizes)

                        if left_intron_size == 0:
                                ORFExonSizes[-1] += IOExonSizes[0]
                                IOExonSizes = [x for idx, x in enumerate(IOExonSizes) if idx != 0]
                                IOIntrons = IOInternalIntrons
                        else:
                                IOIntrons = [left_intron_size] + IOInternalIntrons
                        if (right_intron_size == 0) and (len(IOExonSizes) > 0):
                                CDSExonSizes[0] += IOExonSizes[-1]
                                IOExonSizes = [x for idx, x in enumerate(IOExonSizes) if idx != (len(IOExonSizes) - 1)]
                        else:
                                IOIntrons.append(right_intron_size)
                        AllExons = ORFExonSizes + IOExonSizes + CDSExonSizes
                        AllIntrons = ORFIntronSizes + IOIntrons + CDSIntronSizes
                        
                checkseq = alt_seq + str(interorf)
                ext_checkseq = checkseq + CDSSeq

                if len(checkseq)%3 > 0:
                        ss = Seq(checkseq[0:-(len(checkseq)%3)])
                else:
                        ss = Seq(checkseq)

                if len(ext_checkseq)%3 > 0:
                        ess = Seq(ext_checkseq[0:-(len(ext_checkseq)%3)])
                else:
                        ess = Seq(ext_checkseq)

                sstr = ss.translate()
                esstr = ess.translate()

                if '*' in sstr:
                        main_cds_effect = 'main_CDS_unaffected'                             
                else:
                        #if len(checkseq)%3==0:
                        #print(checkseq)
                        #print(ss)
                        #print(sstr)
                        #print(ext_checkseq)
                        #print(ess)
                        #print(esstr)
                        #print(-(len(ext_checkseq)%3))
                        if (len(ext_checkseq)%3 == 0) and not ('*' in esstr[:-1]):
                                main_cds_effect = 'N-terminal_extension'
                        else:
                                main_cds_effect = 'out-of-frame_overlap'
                
                if '*' not in esstr:
                        ExportExons = AllExons
                        ExportExons[-1] += 3
                        ExportIntrons = AllIntrons
                else:
                        where_stop = [i for i, x in enumerate(esstr) if x == '*'][0]
                        stop_codon_start = (where_stop * 3) + signed_diff
                        ExportExons = []
                        ExportLength = 0
                        ExportIntrons = []
                        for i, k in enumerate(AllExons):
                                if i >= 1:
                                        ExportIntrons.append(AllIntrons[i - 1])
                                if stop_codon_start < (ExportLength + k):
                                        ExportExons.append(stop_codon_start - ExportLength + 2)
                                        break
                                ExportExons.append(k)
                                ExportLength += k

                #print(f"Of a total of {len(AllExons)} remained {len(ExportExons)} due to {conseq} in {uorf_name} on strand {x.loc['strand']} with CDS effect {main_cds_effect}")
                
                AllStarts = np.cumsum([0] + [x + y for x, y in zip(ExportExons[:-1], ExportIntrons)])
                MasterLength = AllStarts[-1] + ExportExons[-1]

                export_bed_line = [x.loc['#CHROM']]

                if x.loc['strand'] == '+':
                        export_bed_line.append(x.loc['orf_start'])
                        export_bed_line.append(x.loc['orf_start'] + MasterLength)
                        ExportStarts = AllStarts
                else:
                        export_bed_line.append(x.loc['orf_end'] - MasterLength)
                        export_bed_line.append(x.loc['orf_end'])
                        ExportStarts = list(reversed([MasterLength - x - y for x, y in zip(ExportExons, AllStarts)]))
                        ExportExons = list(reversed(ExportExons))

                export_name = f"{x.loc['#CHROM']}-{x.loc['POS']}-{x.loc['REF']}>{x.loc['ALT']}|{conseq}|{main_cds_effect}|{uorf_name}"
                export_bed_line.append(export_name)
                export_bed_line.append(0)
                export_bed_line.append(x.loc['strand'])
                export_bed_line.append(export_bed_line[1])
                export_bed_line.append(export_bed_line[1])

                # Color annotation
                if main_cds_effect == 'N-terminal_extension':
                        export_bed_line.append('255,125,0')
                elif main_cds_effect == 'out-of-frame_overlap':
                        export_bed_line.append('255,0,0')
                elif main_cds_effect == 'main_CDS_unaffected':
                        export_bed_line.append('125,125,125')
                else:
                        export_bed_line.append('0,0,0')
                
                export_bed_line.append(len(ExportStarts))
                export_bed_line.append(','.join([str(x) for x in ExportExons]) + ',')
                export_bed_line.append(','.join([str(x) for x in ExportStarts]) + ',')

                export_bed_str = '\t'.join([str(x) for x in export_bed_line])
                #if io_name == 'chr11:31806911-31810840(-)/NM_000280':
                #                 print(export_bed_str)
                if x.loc['codon_type'] == 'ATG' or not atg_only:
                        output_bed.write(f'{export_bed_str}\n')

        # return two columns: symbolic designation of variants and consequence group
        return pd.Series([symb, conseq, main_cds_effect], index=['symbol', 'consequence', 'main_cds_effect'])


if __name__ == '__main__':

        @click.command()
        @click.option('--input_vcf', '-i', required=True)
        @click.option('--bed', '-b', required=True)
        @click.option('--fasta', '-f', required=True)
        @click.option('--gtf', '-g')
        @click.option('--output_tsv', '-ot')
        @click.option('--output_bed', '-ob', required=True)
        @click.option('--output_vcf', '-ov', required=True)
        @click.option('--atg_only', '-atg', is_flag=True, default=False)
        def main(input_vcf, bed, fasta, gtf, output_tsv, output_bed, output_vcf, atg_only) -> None:

                # HARDCODED: annotations from 4 column of bed file
                bed_4col_info = '|utid|overlapping_type|codon_type'
                bed_4col_info_cols = bed_4col_info.split('|')[1:]

                if os.path.isfile(output_bed):
                        print('WARN Output BED file exists, overwriting...')
                out_bed_handle = open(output_bed, 'w')
                out_bed_handle.write('track name="affected uORFs" description=" affected uORFs " itemRgb="On"\n')

                # get VCF-header lines
                tmp_h = NamedTemporaryFile()
                sp.run(f"head -5000 {input_vcf} | grep -P '^##' > {tmp_h.name}", shell=True)
                with open(tmp_h.name) as f:
                        h = f.read()
                # processing
                df = process(input_vcf, bed, fasta, bed_4col_info_cols, gtf, h, out_bed_handle, atg_only)
                if atg_only:
                        df = df.loc[df['codon_type'] == 'ATG']

                df.drop_duplicates(subset=['#CHROM', 'POS', 'REF', 'ALT',
                                        'orf_start', 'orf_end'], inplace=True)

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
                        f'{x["main_cds_effect"]}|'
                        f'{x["in_known_CDS"]}|'
                        f'{x["in_known_ORF"]}' for x in df.to_dict(orient='records')]

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
                df['INFO'] = df['INFO'].astype(str) + ';' + df['INFO_new'].astype(str)
                df = df.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]

                # remove output file before appending new data
                out_vcf = Path(output_vcf)
                if out_vcf.exists():
                        out_vcf.unlink()

                # add header ##INFO uORF_annotator (uBERT) line
                h += \
                f'##INFO=<ID=uBERT,Number=.,Type=String,Description="Consequence uORF_annotator from uBERT. '
                f'Format: ORF_START|ORF_END|ORF_SYMB|ORF_CONSEQ|main_cds_effect|in_known_CDS|in_known_ORF{bed_4col_info}">'
                # write VCF-header and VCF-body in output file
                with open(output_vcf, 'a') as w:
                        w.write(h)
                        w.write('\n')
                        df.to_csv(w, sep='\t', index=None)
                
                out_bed_handle.close()
        # run analysis
        main()