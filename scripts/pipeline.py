from data_processor import RawDataProcessor, uORFDataProcessor, InterORFDataProcessor


class DataProcessor:
    def __init__(self, gtf_file_path: str, vcf_file_path: str, bed_file_path: str, fasta_file_path: str):
        self.raw_data_processor: RawDataProcessor = RawDataProcessor()
        self.uorf_processor: uORFDataProcessor = uORFDataProcessor()
        self.interorf_processor: InterORFDataProcessor = InterORFDataProcessor()

        self.gtf: str = gtf_file_path
        self.vcf: str = vcf_file_path
        self.bed: str = bed_file_path
        self.fasta: str = fasta_file_path

        self.bed_4col_info: str = '|utid|overlapping_type|dominance_type|codon_type'
        self.bed_4col_info_cols: list[str] = self.bed_4col_info.split('|')[1:]

    def process_data(self) -> None:
        self._process_raw_data(self.vcf, self.bed, self.gtf, self.fasta)
        self._process_uorf_data(self.raw_data_processor, self.bed_4col_info_cols)
        self._process_interorf_data(self.raw_data_processor, self.uorf_processor)

    def _process_raw_data(self, vcf_file_path: str, bed_file_path: str, gtf_file_path: str, fasta_file_path: str) -> None:
        self.raw_data_processor.process_data(vcf_file_path, bed_file_path, gtf_file_path, fasta_file_path)

    def _process_uorf_data(self, raw_data_processor: RawDataProcessor, bed_4col_info_cols: list[str]) -> None:
        self.uorf_processor.process_data(bed_4col_info_cols, raw_data_processor)

    def _process_interorf_data(self, raw_data_processor: RawDataProcessor, uorf_processor: uORFDataProcessor) -> None:
        self.interorf_processor.process_data(raw_data_processor, uorf_processor)
