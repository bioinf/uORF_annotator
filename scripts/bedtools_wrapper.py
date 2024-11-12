from pybedtools import BedTool


class Bedtools:
	@staticmethod
	def intersect(file1, file2):
		bed_file1 = BedTool(file1)
		bed_file2 = BedTool(file2)
		intersected_bed = bed_file1.intersect(bed_file2, wo=True)
		intersection_df = intersected_bed.to_dataframe(header=None)
		intersection_df.columns = range(len(intersection_df.columns))
		return intersection_df
	
	@staticmethod
	def getfasta(fasta_file, bed_file, strand=True):
		fasta = BedTool(fasta_file)
		bed = BedTool(bed_file)
		getfasta = bed.sequence(fi=fasta, s=strand)
		return getfasta.seqfn
