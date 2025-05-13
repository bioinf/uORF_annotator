import logging
from typing import Optional, Dict, List, Set, Tuple
import pandas as pd
import pysam
import hashlib
import json

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact
from transcript_sequence import TranscriptSequence

class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile, bed_output: str = None, debug_mode: bool = False):
        """Initialize the variant processor with coordinate converter and reference sequence."""
        self.converter = converter
        self.fasta = fasta
        self.bed_output = bed_output
        self.debug_mode = debug_mode
        self.annotator = None
        self.debug_info = {}

    def process_variant(self, row: pd.Series) -> List[Dict]:
        """
        Process a variant to determine its effect on uORF and main CDS.
        Returns a list of results for all variant-uORF pairs, avoiding redundancy
        between regular and extended transcripts.
        """
        try:
            # Extract variant information
            chrom = row['col0']
            vcf_pos = int(row['col1'])
            
            # Try to extract the RSID (careful with different VCF formats)
            try:
                rsid = str(row.get('col2', '.'))
                if rsid == 'nan' or rsid == 'None':
                    rsid = '.'
            except:
                rsid = '.'
                
            ref_allele = row['col3']
            alt_allele = row['col4']
            
            # Check if we have combined BED entries from the pipeline
            if 'all_bed_entries' in row:
                bed_entries = row['all_bed_entries']
                logging.info(f"Processing variant {chrom}:{vcf_pos} {ref_allele}>{alt_allele} with {len(bed_entries)} BED entries")
            else:
                # Use the single BED entry in the row
                bed_full_name = row['col11']
                bed_entries = [{
                    'start': int(row['col9']),
                    'end': int(row['col10']),
                    'name': bed_full_name,
                    'strand': row['col13']
                }]
                logging.info(f"Processing variant {chrom}:{vcf_pos} {ref_allele}>{alt_allele} with single BED entry: {bed_full_name}")
            
            # Reset tracking for this variant
            self._processed_uorfs = {}
            
            # Process each BED entry
            results = []
            for bed_entry in bed_entries:
                bed_full_name = bed_entry['name']
                bed_start = bed_entry['start']
                bed_end = bed_entry['end']
                
                # Extract transcript ID
                original_transcript_id = self._extract_transcript_id(bed_full_name)
                
                # Find all matching transcripts for this BED entry
                matching_transcripts = []
                for tid, transcript_obj in self.converter.transcripts.items():
                    # Check if this transcript corresponds to this BED entry
                    if (tid == original_transcript_id or tid.startswith(f"{original_transcript_id}_uorf_")) and \
                       transcript_obj.uorf_start_genomic == bed_start + 1 and \
                       transcript_obj.uorf_end_genomic == bed_end:
                        matching_transcripts.append(tid)
                
                if not matching_transcripts:
                    logging.warning(f"No matching transcripts found for BED entry: {bed_full_name}")
                    continue
                
                if self.debug_mode:
                    logging.debug(f"Found {len(matching_transcripts)} matching transcripts for {bed_full_name}: {matching_transcripts}")
                
                # Separate regular and extended transcripts
                regular_transcripts = []
                extended_transcripts = []
                
                for tid in matching_transcripts:
                    transcript_obj = self.converter.transcripts[tid]
                    was_extended = getattr(transcript_obj, 'was_extended', False)
                    
                    if was_extended:
                        extended_transcripts.append(tid)
                    else:
                        regular_transcripts.append(tid)
                
                # Track positions already processed
                processed_positions = set()
                
                # First try to process regular transcripts
                for transcript_id in regular_transcripts:
                    result = self._process_single_transcript(
                        transcript_id=transcript_id,
                        chrom=chrom,
                        vcf_pos=vcf_pos,
                        ref_allele=ref_allele,
                        alt_allele=alt_allele,
                        rsid=rsid,
                        bed_full_name=bed_full_name,
                        processed_positions=processed_positions,
                        is_extended=False
                    )
                    if result:
                        # Add result and update processed_positions
                        if 'Transcript_Position' in result:
                            processed_positions.add(result['Transcript_Position'])
                        
                        # Check for duplicates before adding
                        result_id = self._create_result_identifier(result)
                        if not any(result_id == self._create_result_identifier(r) for r in results):
                            results.append(result)
                            if self.debug_mode:
                                logging.debug(f"Added result for transcript-uORF: {transcript_id}")
                        else:
                            if self.debug_mode:
                                logging.debug(f"Skipping duplicate result for transcript-uORF: {transcript_id}")
                
                # If no results from regular transcripts, try extended transcripts
                if not results and extended_transcripts:
                    if self.debug_mode:
                        logging.debug("No results from regular transcripts, trying extended transcripts")
                    
                    for transcript_id in extended_transcripts:
                        result = self._process_single_transcript(
                            transcript_id=transcript_id,
                            chrom=chrom,
                            vcf_pos=vcf_pos,
                            ref_allele=ref_allele,
                            alt_allele=alt_allele,
                            rsid=rsid,
                            bed_full_name=bed_full_name,
                            processed_positions=processed_positions,
                            is_extended=True
                        )
                        if result:
                            # Check for duplicates before adding
                            result_id = self._create_result_identifier(result)
                            if not any(result_id == self._create_result_identifier(r) for r in results):
                                results.append(result)
                                if self.debug_mode:
                                    logging.debug(f"Added result for extended transcript-uORF: {transcript_id}")
                            else:
                                if self.debug_mode:
                                    logging.debug(f"Skipping duplicate result for extended transcript-uORF: {transcript_id}")
                
            # Return all results
            return results
                    
        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}", exc_info=True)
            return []

    def _process_single_transcript(self, transcript_id, chrom, vcf_pos, ref_allele, alt_allele, 
                                  rsid, bed_full_name, processed_positions, is_extended):
        """
        Process a single transcript to determine variant effect.
        This consolidated method handles both regular and extended transcripts.
        
        Args:
            transcript_id: The transcript ID to process
            chrom: Chromosome name
            vcf_pos: Variant position in VCF (genomic) coordinates
            ref_allele: Reference allele
            alt_allele: Alternative allele
            rsid: Variant ID from VCF
            bed_full_name: Full name from BED file
            processed_positions: Set of already processed transcript positions
            is_extended: Whether this is an extended transcript
            
        Returns:
            Result dictionary or None if processing fails
        """
        transcript_obj = self.converter.transcripts[transcript_id]
        
        # Create a unique identifier for this uORF based on genomic coordinates
        uorf_key = self._create_uorf_key(transcript_obj)
        
        # Skip if we've already processed this uORF
        if uorf_key in self._processed_uorfs:
            if self.debug_mode:
                logging.debug(f"Skipping already processed uORF: {uorf_key}")
            return None
            
        # Mark this uORF as processed
        self._processed_uorfs[uorf_key] = True
        
        # Log debug information
        if self.debug_mode and is_extended:
            self._log_debug_info(transcript_id, transcript_obj, vcf_pos)
        
        # Check if the variant is actually within the uORF genomic coordinates
        if (transcript_obj.uorf_start_genomic is not None and 
            transcript_obj.uorf_end_genomic is not None):
            
            # For positive strand, check if position is within uORF boundaries
            if transcript_obj.strand == '+':
                is_within_uorf = (vcf_pos >= transcript_obj.uorf_start_genomic and 
                              vcf_pos <= transcript_obj.uorf_end_genomic)
            # For negative strand, need to handle possible swapped coordinates
            else:
                uorf_min = min(transcript_obj.uorf_start_genomic, transcript_obj.uorf_end_genomic)
                uorf_max = max(transcript_obj.uorf_start_genomic, transcript_obj.uorf_end_genomic)
                is_within_uorf = (vcf_pos >= uorf_min and vcf_pos <= uorf_max)
            
            if not is_within_uorf:
                if self.debug_mode:
                    logging.debug(f"Position {vcf_pos} is outside uORF genomic coordinates "
                               f"({transcript_obj.uorf_start_genomic}-{transcript_obj.uorf_end_genomic}) "
                               f"for {'extended' if is_extended else ''} transcript {transcript_id}, skipping")
                return None  # Skip to next transcript
        
        # Get transcript position for the variant
        variant_coords = transcript_obj.get_coordinates(vcf_pos)
        if not variant_coords:
            if self.debug_mode:
                logging.debug(f"Position {vcf_pos} not in transcript coordinates for {transcript_id}, skipping")
            return None  # Skip to next transcript
        
        # Check if this position was already processed in a regular transcript
        if variant_coords.transcript in processed_positions:
            if self.debug_mode:
                logging.debug(f"Position {variant_coords.transcript} already processed in regular transcript, skipping")
            return None
            
        # Check if uORF coordinates are still missing after potential extension
        if transcript_obj.uorf_start is None or transcript_obj.uorf_end is None:
            if self.debug_mode:
                logging.debug(f"Missing uORF transcript coordinates for {transcript_id}, skipping")
            return None  # Skip to next transcript

        # Handle variant alleles based on strand
        # For negative strand, we need to reverse-complement the ref and alt alleles
        variant_ref = ref_allele
        variant_alt = alt_allele
        
        if transcript_obj.strand == '-':
            # Use the TranscriptSequence's reverse complement method
            variant_ref = TranscriptSequence._reverse_complement(ref_allele)
            variant_alt = TranscriptSequence._reverse_complement(alt_allele)
            if self.debug_mode:
                logging.debug(f"Converted alleles for negative strand: {ref_allele}>{alt_allele} to {variant_ref}>{variant_alt}")
            
        # Create TranscriptSequence object
        transcript_seq = TranscriptSequence(transcript_obj, self.fasta, chrom, debug_mode=self.debug_mode)
        if not transcript_seq.sequence:
            if self.debug_mode:
                logging.debug(f"Could not extract transcript sequence for {transcript_id}, skipping")
            return None
            
        # Check if uORF region was extracted successfully
        if not transcript_seq.uorf_region:
            if self.debug_mode:
                logging.debug(f"Failed to extract uORF region for {transcript_id}, skipping")
            return None
        
        # Get start codon information
        start_codon, is_canonical_atg = transcript_seq.get_start_codon()
        start_codon_type = "ATG" if is_canonical_atg else "NON-ATG"
        
        # Update transcript object with verified start codon type
        transcript_obj.start_codon_type = start_codon_type
            
        # Initialize annotator with TranscriptSequence
        self.annotator = VariantAnnotator(transcript_seq, transcript_obj, bed_file_path=self.bed_output, debug_mode=self.debug_mode)

        # Get pre-determined overlap status from the Transcript object
        overlaps_maincds = getattr(transcript_obj, 'overlaps_maincds', False)
        if self.debug_mode:
            logging.debug(f"Using pre-determined overlap status for {transcript_id}: {overlaps_maincds}")

        # Get codon change
        codon_change = self.annotator.get_codon_change({
            'position': variant_coords.transcript,
            'ref_allele': variant_ref,
            'alt_allele': variant_alt
        })

        # Log detailed information about codon processing
        if self.debug_mode:
            logging.debug(f"Codon processing details for {transcript_id}:")
            logging.debug(f"  Original alleles: {ref_allele}>{alt_allele}")
            if transcript_obj.strand == '-':
                logging.debug(f"  Strand: negative, using reverse complement")
                logging.debug(f"  Reverse complemented alleles: {variant_ref}>{variant_alt}")
            else:
                logging.debug(f"  Strand: positive, using original alleles")
                logging.debug(f"  Processed alleles: {variant_ref}>{variant_alt}")
            logging.debug(f"  Computed codon change: {codon_change}")

        # Prepare variant data with all information
        variant_data = {
            'chromosome': chrom,
            'rsid': rsid,
            'full_uorf_name': bed_full_name,
            'position': variant_coords.transcript,
            'position_genomic': vcf_pos,
            'ref_allele': variant_ref,
            'alt_allele': variant_alt,
            'strand': transcript_obj.strand,
            'uorf_start': transcript_obj.uorf_start,
            'uorf_end': transcript_obj.uorf_end,
            'maincds_start': transcript_obj.mainorf_start,
            'maincds_end': transcript_obj.mainorf_end,
            'uorf_start_genomic': transcript_obj.uorf_start_genomic,
            'uorf_end_genomic': transcript_obj.uorf_end_genomic,
            'maincds_start_genomic': transcript_obj.mainorf_start_genomic,
            'maincds_end_genomic': transcript_obj.mainorf_end_genomic,
            'codon_change': codon_change,
            'overlaps_maincds': overlaps_maincds,
            'transcript_extended': is_extended,  # This should reflect extended status
            'start_codon_type': start_codon_type,  # Add start codon type information
            'start_codon': start_codon  # Add actual start codon sequence
        }

        # Get consequences and impacts
        uorf_consequence = self.annotator.get_consequence(variant_data)
        maincds_impact = None
        if uorf_consequence:
            maincds_impact = self.annotator.predict_impact(variant_data, uorf_consequence)
        
        # Get corrected uORF coordinates that reflect biological reality
        uorf_start_transcript, uorf_end_transcript = self._get_corrected_uorf_coordinates(
            transcript_obj, transcript_seq
        )

        # Get corrected mainCDS coordinates
        maincds_start_transcript = transcript_obj.mainorf_start
        maincds_end_transcript = transcript_obj.mainorf_end

        # Final check and cleanup - ensure we have valid numerical values
        uorf_start_transcript = uorf_start_transcript if uorf_start_transcript is not None else -1
        uorf_end_transcript = uorf_end_transcript if uorf_end_transcript is not None else -1
        maincds_start_transcript = maincds_start_transcript if maincds_start_transcript is not None else -1
        maincds_end_transcript = maincds_end_transcript if maincds_end_transcript is not None else -1

        # Extract the original transcript ID if this is a derived uORF transcript
        display_transcript_id = self._extract_transcript_id(bed_full_name)
        
        result = {
            'Chromosome': chrom,
            'Original_Genome_Position': vcf_pos,
            'RSID': rsid,
            'BED_Name': bed_full_name,
            'Transcript_Position': variant_coords.transcript,
            'Ref_Allele': ref_allele,  # Save original alleles for output
            'Alt_Allele': alt_allele,
            'Strand': transcript_obj.strand,
            'Transcript_ID': display_transcript_id,
            'uORF_Start_Genomic': transcript_obj.uorf_start_genomic,
            'uORF_End_Genomic': transcript_obj.uorf_end_genomic,
            'uORF_Start_Transcript': uorf_start_transcript,
            'uORF_End_Transcript': uorf_end_transcript,
            'mainCDS_Start_Genomic': transcript_obj.mainorf_start_genomic,
            'mainCDS_End_Genomic': transcript_obj.mainorf_end_genomic,
            'mainCDS_Start_Transcript': maincds_start_transcript,
            'mainCDS_End_Transcript': maincds_end_transcript,
            'Codon_Change': codon_change,
            'uORF_Consequence': uorf_consequence.value if uorf_consequence else 'None',
            'uORF_mainCDS_Overlap': 'overlapping' if overlaps_maincds else 'non_overlapping',
            'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None',
            'Start_Codon_Type': start_codon_type,  # Add start codon type to the result
            'Start_Codon': start_codon  # Add actual start codon sequence to the result
        }
        
        # Add additional information if transcript was extended
        if is_extended:
            result['Transcript_Extended'] = 'Yes'
        
        return result
    
    def _create_result_identifier(self, result: Dict) -> str:
        """Create a unique identifier for a result to detect duplicates."""
        # The key fields that make a result unique
        key_fields = [
            result['Chromosome'],
            str(result['Original_Genome_Position']),
            result['Ref_Allele'],
            result['Alt_Allele'],
            result['Transcript_ID'],
            str(result['uORF_Start_Genomic']),
            str(result['uORF_End_Genomic']),
            result['Codon_Change'],
            result['uORF_Consequence'],
            result['mainCDS_Impact']
        ]
        
        # Join fields with a delimiter
        return '|'.join(key_fields)
    
    def _create_uorf_key(self, transcript_obj) -> str:
        """
        Create a unique key for a uORF based on its genomic coordinates.
        This helps identify and prevent processing duplicate uORFs.
        """
        if transcript_obj.uorf_start_genomic is None or transcript_obj.uorf_end_genomic is None:
            return f"{transcript_obj.transcript_id}"
            
        return f"{transcript_obj.chromosome}:{transcript_obj.uorf_start_genomic}-{transcript_obj.uorf_end_genomic}"
    
    def _get_corrected_uorf_coordinates(self, transcript_obj, transcript_seq=None):
            """
            Get corrected uORF coordinates for output.
            
            Args:
                transcript_obj: Transcript object containing coordinate information
                transcript_seq: Optional TranscriptSequence object with stored original coordinates
                
            Returns:
                Tuple of (uorf_start_transcript, uorf_end_transcript)
            """
            # Use coordinates from Transcript directly
            uorf_start = transcript_obj.uorf_start
            uorf_end = transcript_obj.uorf_end
            
            # Validate coordinates
            if uorf_start is None or uorf_end is None:
                logging.warning(f"Missing uORF transcript coordinates for {transcript_obj.transcript_id}")
                # Use genomic coordinates as fallbacks if needed
                if uorf_start is None and transcript_obj.uorf_start_genomic is not None:
                    uorf_start = transcript_obj.get_transcript_position(transcript_obj.uorf_start_genomic)
                    if self.debug_mode:
                        logging.debug(f"Using fallback for uORF start: {uorf_start}")
                if uorf_end is None and transcript_obj.uorf_end_genomic is not None:
                    uorf_end = transcript_obj.get_transcript_position(transcript_obj.uorf_end_genomic)
                    if self.debug_mode:
                        logging.debug(f"Using fallback for uORF end: {uorf_end}")
            
            return uorf_start, uorf_end

    def _log_debug_info(self, transcript_id, transcript_obj, vcf_pos):
        """Log debug information about transcript coordinates."""
        debug_info = {
            'transcript_id': transcript_id,
            'strand': transcript_obj.strand,
            'uorf_start_g': transcript_obj.uorf_start_genomic,
            'uorf_end_g': transcript_obj.uorf_end_genomic,
            'uorf_start_t': transcript_obj.uorf_start,
            'uorf_end_t': transcript_obj.uorf_end,
            'variant_pos_g': vcf_pos,
            'was_extended': getattr(transcript_obj, 'was_extended', False),
            'overlaps_maincds': getattr(transcript_obj, 'overlaps_maincds', False),
            'start_codon_type': getattr(transcript_obj, 'start_codon_type', 'Unknown')
        }
        
        # Store debug info
        self.debug_info[transcript_id] = debug_info
        
        # Check for potential issues
        if transcript_obj.uorf_start is None and transcript_obj.uorf_start_genomic is not None:
            logging.warning(f"uORF genomic coordinates exist but transcript coordinates are None for {transcript_id}")
            
        # Show transcript boundaries
            if transcript_obj.genome_to_transcript:
                min_pos = min(transcript_obj.genome_to_transcript.keys())
                max_pos = max(transcript_obj.genome_to_transcript.keys())
                
                # Check if uORF is outside transcript bounds
                if transcript_obj.uorf_start_genomic < min_pos or transcript_obj.uorf_start_genomic > max_pos:
                    logging.warning(f"uORF start ({transcript_obj.uorf_start_genomic}) is outside transcript bounds ({min_pos}-{max_pos})")
                if transcript_obj.uorf_end_genomic < min_pos or transcript_obj.uorf_end_genomic > max_pos:
                    logging.warning(f"uORF end ({transcript_obj.uorf_end_genomic}) is outside transcript bounds ({min_pos}-{max_pos})")
        
        if transcript_obj.uorf_start is not None and transcript_obj.uorf_end is not None:
            if transcript_obj.strand == '+' and transcript_obj.uorf_start > transcript_obj.uorf_end:
                logging.warning(f"Suspicious positive strand coordinates: uORF_start > uORF_end for {transcript_id}")
            elif transcript_obj.strand == '-' and transcript_obj.uorf_start < transcript_obj.uorf_end:
                logging.warning(f"Suspicious negative strand coordinates: uORF_start < uORF_end for {transcript_id}")
            
            # Check for invalid coordinates
            if transcript_obj.uorf_start < 1 or transcript_obj.uorf_end < 1:
                logging.warning(f"Invalid transcript coordinates: less than 1 for {transcript_id}")

    @staticmethod
    def _extract_transcript_id(bed_name: str) -> str:
        """Extract transcript ID from BED name field."""
        return bed_name.split('|')[0].split('.')[0]