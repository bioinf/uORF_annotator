import logging
from typing import Optional, Dict
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact


class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile):
        self.converter = converter
        self.fasta = fasta
        self.annotator = None

    def _get_sequence_context(self, chrom: str, start: int, end: int) -> str:
        try:
            return self.fasta.fetch(chrom, start, end).upper()
        except Exception as e:
            logging.error(f"Error fetching sequence for {chrom}:{start}-{end}: {e}")
            return ""

    def extract_uorf_sequence(self, transcript_obj, chrom: str) -> str:
        if transcript_obj.uorf_start_genomic is None or transcript_obj.uorf_end_genomic is None:
            logging.error("Missing uORF genomic coordinates")
            return ""
            
        uorf_start = transcript_obj.uorf_start_genomic
        uorf_end = transcript_obj.uorf_end_genomic
        uorf_seq = ""

        uorf_exons_to_process = (
            transcript_obj.exons if transcript_obj.strand == '+' 
            else reversed(transcript_obj.exons)
        )

        for i, exon in enumerate(uorf_exons_to_process):
            if exon.genome_end < uorf_start or exon.genome_start > uorf_end:
                continue
            
            start = max(exon.genome_start, uorf_start)
            end = min(exon.genome_end, uorf_end)
            
            try:
                if transcript_obj.strand == '+':
                    exon_seq = self._get_sequence_context(
                        chrom, 
                        start if i == 0 else start - 1,
                        end
                    )
                else:
                    exon_seq = self._get_sequence_context(
                        chrom, 
                        start if i == 0 else start - 1,
                        end if i == 0 else end - 1
                    )
                    exon_seq = self._reverse_complement(exon_seq)
                
                uorf_seq += exon_seq
            except Exception as e:
                logging.error(f"Error extracting sequence for exon: {str(e)}")
                continue

        return uorf_seq

    def process_variant(self, row: pd.Series) -> Optional[Dict]:
            try:
                # Clear previous annotator
                self.annotator = None

                # Extract variant information
                chrom = row['col0']
                vcf_pos = int(row['col1'])
                ref_allele = row['col3']
                alt_allele = row['col4']
                uorf_start_genomic = int(row['col9'])
                uorf_end_genomic = int(row['col10'])
                strand = row['col13']
                
                # Get transcript information
                transcript_id = self._extract_transcript_id(row['col11'])
                matching_transcripts = [
                    tid for tid in self.converter.transcripts.keys() 
                    if tid == transcript_id
                ]

                if not matching_transcripts:
                    logging.warning(f"Transcript {transcript_id} not found")
                    return None

                matched_transcript = matching_transcripts[0]
                transcript_obj = self.converter.transcripts[matched_transcript]

                # Get transcript position
                variant_coords = transcript_obj.get_coordinates(vcf_pos)
                if not variant_coords:
                    logging.warning(f"Could not convert position {vcf_pos} to transcript coordinates")
                    return None

                if transcript_obj.uorf_start is None or transcript_obj.uorf_end is None:
                    logging.warning(f"Missing uORF transcript coordinates for {transcript_id}")
                    return None

                # Extract sequence and initialize annotator
                uorf_seq = self.extract_uorf_sequence(transcript_obj, chrom)
                if not uorf_seq:
                    logging.warning(f"Could not extract uORF sequence for {transcript_id}")
                    return None
                    
                self.annotator = VariantAnnotator(uorf_seq, self.fasta)
                
                # Prepare variant data with validation
                variant_data = {
                    'position': variant_coords.transcript,
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele,
                    'uorf_start': transcript_obj.uorf_start,
                    'uorf_end': transcript_obj.uorf_end,
                    'strand': strand,
                    'maincds_start': transcript_obj.mainorf_start,
                    'maincds_end': transcript_obj.mainorf_end,
                    'chromosome': chrom,
                    'genomic_start': uorf_start_genomic,
                    'genomic_end': uorf_end_genomic
                }

                # Get consequences and impacts
                uorf_consequence = self.annotator.get_consequence(variant_data)
                maincds_impact = None
                if uorf_consequence:
                    maincds_impact = self.annotator.predict_impact(variant_data, uorf_consequence)
                
                # Get codon change
                codon_change = self.annotator.get_codon_change(variant_data)
                
                # Check overlap
                overlaps_maincds = (
                    self.annotator.does_overlap_maincds(
                        transcript_obj.uorf_start,
                        transcript_obj.uorf_end,
                        transcript_obj.mainorf_start,
                        transcript_obj.mainorf_end,
                        strand
                    )
                    if (transcript_obj.uorf_start is not None and 
                        transcript_obj.uorf_end is not None and 
                        transcript_obj.mainorf_start is not None and 
                        transcript_obj.mainorf_end is not None)
                    else False
                )
                
                return {
                    'Chromosome': chrom,
                    'Original_Genome_Position': vcf_pos,
                    'Transcript_Position': variant_coords.transcript,
                    'Ref_Allele': ref_allele,
                    'Alt_Allele': alt_allele,
                    'Strand': strand,
                    'Transcript_ID': matched_transcript,
                    'uORF_Start_Genomic': transcript_obj.uorf_start_genomic,
                    'uORF_End_Genomic': transcript_obj.uorf_end_genomic,
                    'uORF_Start_Transcript': transcript_obj.uorf_start,
                    'uORF_End_Transcript': transcript_obj.uorf_end,
                    'mainCDS_Start_Genomic': transcript_obj.mainorf_start_genomic,
                    'mainCDS_End_Genomic': transcript_obj.mainorf_end_genomic,
                    'mainCDS_Start_Transcript': transcript_obj.mainorf_start,
                    'mainCDS_End_Transcript': transcript_obj.mainorf_end,
                    'Codon_Change': codon_change,
                    'uORF_Consequence': uorf_consequence.value if uorf_consequence else 'None',
                    'uORF_mainCDS_Overlap': 'overlapping' if overlaps_maincds else 'non_overlapping',
                    'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None'
                }
                
            except Exception as e:
                logging.error(f"Error processing variant: {str(e)}")
                return None

    def _reverse_complement(self, sequence: str) -> str:
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base.upper(), base) for base in reversed(sequence))

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        return bed_name.split('|')[0].split('.')[0]