from pybedtools import BedTool
import pandas as pd


class Bedtools:
	@staticmethod
	def intersect(file1: str, file2: str) -> pd.DataFrame:
		bed_file1 = BedTool(file1)
		bed_file2 = BedTool(file2)
		intersected_bed = bed_file1.intersect(bed_file2, wo=True)
		intersection_df = intersected_bed.to_dataframe(header=None)
		intersection_df.columns = range(len(intersection_df.columns))
		return intersection_df
	
	@staticmethod
	def getfasta(fasta_file: str, bed_file: str, strand: bool = True) -> str:
		fasta = BedTool(fasta_file)
		bed = BedTool(bed_file)
		getfasta = bed.sequence(fi=fasta, s=strand)
		return getfasta.seqfn
