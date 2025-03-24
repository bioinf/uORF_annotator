import logging
from typing import Optional, Dict, List
import pandas as pd
import pysam

from converters import CoordinateConverter
from annotator import VariantAnnotator, UORFConsequence, MainCDSImpact
from transcript_sequence import TranscriptSequence

class VariantProcessor:
    def __init__(self, converter: CoordinateConverter, fasta: pysam.FastaFile):
        """Initialize the variant processor with coordinate converter and reference sequence."""
        self.converter = converter
        self.fasta = fasta
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
            ref_allele = row['col3']
            alt_allele = row['col4']
            
            logging.debug(f"Processing variant at {chrom}:{vcf_pos} {ref_allele}>{alt_allele}")
            
            # Get transcript information
            original_transcript_id = self._extract_transcript_id(row['col11'])
            
            # Find all matching transcripts, including those with multiple uORFs
            matching_transcripts = [
                tid for tid in self.converter.transcripts.keys() 
                if tid == original_transcript_id or tid.startswith(f"{original_transcript_id}_uorf_")
            ]
            
            if not matching_transcripts:
                logging.warning(f"No transcripts found for {original_transcript_id}")
                return []
            
            logging.debug(f"Found {len(matching_transcripts)} matching transcripts: {matching_transcripts}")
            
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
            
            logging.debug(f"Regular transcripts: {regular_transcripts}")
            logging.debug(f"Extended transcripts: {extended_transcripts}")
            
            # Process each matching transcript-uORF pair
            results = []
            processed_positions = set()  # Track positions already processed
            
            # First try to process regular transcripts
            for transcript_id in regular_transcripts:
                logging.debug(f"Processing regular transcript-uORF: {transcript_id}")
                
                transcript_obj = self.converter.transcripts[transcript_id]
                
                # Log debug information
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
                        logging.debug(f"Position {vcf_pos} is outside uORF genomic coordinates "
                                    f"({transcript_obj.uorf_start_genomic}-{transcript_obj.uorf_end_genomic}) "
                                    f"for {transcript_id}, skipping")
                        continue  # Skip to next transcript
                
                # Get transcript position for the variant
                variant_coords = transcript_obj.get_coordinates(vcf_pos)
                if variant_coords:
                    # If we found coordinates in a regular transcript, add to processed positions
                    processed_positions.add(variant_coords.transcript)
                    logging.debug(f"Found position {vcf_pos} in regular transcript {transcript_id}")
                else:
                    logging.debug(f"Position {vcf_pos} not in regular transcript coordinates for {transcript_id}, skipping")
                    continue  # Skip to next transcript
                    
                # Check if uORF coordinates are still missing after potential extension
                if transcript_obj.uorf_start is None or transcript_obj.uorf_end is None:
                    logging.debug(f"Missing uORF transcript coordinates for {transcript_id}, skipping")
                    continue  # Skip to next transcript
                
                # Handle variant alleles based on strand
                # For negative strand, we need to reverse-complement the ref and alt alleles
                # because we're working with sequences in transcript orientation (5' to 3')
                variant_ref = ref_allele
                variant_alt = alt_allele
                
                if transcript_obj.strand == '-':
                    # Use the TranscriptSequence's reverse complement method since we need
                    # to transform alleles to match the transcript orientation
                    variant_ref = TranscriptSequence._reverse_complement(ref_allele)
                    variant_alt = TranscriptSequence._reverse_complement(alt_allele)
                    logging.debug(f"Converted alleles for negative strand: {ref_allele}>{alt_allele} to {variant_ref}>{variant_alt}")

                # Create TranscriptSequence object
                transcript_seq = TranscriptSequence(transcript_obj, self.fasta, chrom)
                if not transcript_seq.sequence:
                    logging.debug(f"Could not extract transcript sequence for {transcript_id}, skipping")
                    continue
                    
                # Check if uORF region was extracted successfully
                if not transcript_seq.uorf_region:
                    logging.debug(f"Failed to extract uORF region for {transcript_id}, skipping")
                    continue
                    
                # Initialize annotator with TranscriptSequence
                self.annotator = VariantAnnotator(transcript_seq)

                # Get pre-determined overlap status from the Transcript object
                # This was determined at transcript creation time using raw coordinates
                overlaps_maincds = getattr(transcript_obj, 'overlaps_maincds', False)
                logging.debug(f"Using pre-determined overlap status for {transcript_id}: {overlaps_maincds}")

                # Get codon change
                codon_change = self.annotator.get_codon_change({
                    'position': variant_coords.transcript,
                    'ref_allele': variant_ref,
                    'alt_allele': variant_alt
                })
                    
                # Prepare variant data with all information
                variant_data = {
                    'position': variant_coords.transcript,
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
                    'transcript_extended': getattr(transcript_obj, 'was_extended', False)
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
                display_transcript_id = original_transcript_id
                    
                result = {
                    'Chromosome': chrom,
                    'Original_Genome_Position': vcf_pos,
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
                    'mainCDS_Impact': maincds_impact.value if maincds_impact else 'None'
                }
                
                # Add additional information if transcript was extended
                was_extended = getattr(transcript_obj, 'was_extended', False)
                if was_extended:
                    result['Transcript_Extended'] = 'Yes'
                    
                results.append(result)
                logging.debug(f"Successfully processed variant for transcript-uORF: {transcript_id}")
                
                # Add this position to processed positions so we don't process it again in extended transcripts
                processed_positions.add(variant_coords.transcript)
                    
            # If no results from regular transcripts, try extended transcripts
            if not results and extended_transcripts:
                logging.debug("No results from regular transcripts, trying extended transcripts")
                
                for transcript_id in extended_transcripts:
                    logging.debug(f"Processing extended transcript-uORF: {transcript_id}")
                    
                    transcript_obj = self.converter.transcripts[transcript_id]
                    
                    # Log debug information
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
                            logging.debug(f"Position {vcf_pos} is outside uORF genomic coordinates "
                                        f"({transcript_obj.uorf_start_genomic}-{transcript_obj.uorf_end_genomic}) "
                                        f"for extended transcript {transcript_id}, skipping")
                            continue  # Skip to next transcript
                    
                    # Get transcript position for the variant
                    variant_coords = transcript_obj.get_coordinates(vcf_pos)
                    if not variant_coords:
                        logging.debug(f"Position {vcf_pos} not in extended transcript coordinates for {transcript_id}, skipping")
                        continue  # Skip to next transcript
                    
                    # Check if this position was already processed in a regular transcript
                    if variant_coords.transcript in processed_positions:
                        logging.debug(f"Position {variant_coords.transcript} already processed in regular transcript, skipping")
                        continue
                        
                    # Check if uORF coordinates are still missing after potential extension
                    if transcript_obj.uorf_start is None or transcript_obj.uorf_end is None:
                        logging.debug(f"Missing uORF transcript coordinates for {transcript_id}, skipping")
                        continue  # Skip to next transcript

                    # Handle variant alleles based on strand
                    # For negative strand, we need to reverse-complement the ref and alt alleles
                    variant_ref = ref_allele
                    variant_alt = alt_allele
                    
                    if transcript_obj.strand == '-':
                        # Use the TranscriptSequence's reverse complement method
                        variant_ref = TranscriptSequence._reverse_complement(ref_allele)
                        variant_alt = TranscriptSequence._reverse_complement(alt_allele)
                        logging.debug(f"Converted alleles for negative strand: {ref_allele}>{alt_allele} to {variant_ref}>{variant_alt}")
                        
                    # Create TranscriptSequence object
                    transcript_seq = TranscriptSequence(transcript_obj, self.fasta, chrom)
                    if not transcript_seq.sequence:
                        logging.debug(f"Could not extract transcript sequence for {transcript_id}, skipping")
                        continue
                        
                    # Check if uORF region was extracted successfully
                    if not transcript_seq.uorf_region:
                        logging.debug(f"Failed to extract uORF region for {transcript_id}, skipping")
                        continue
                        
                    # Initialize annotator with TranscriptSequence
                    self.annotator = VariantAnnotator(transcript_seq)

                    # Get pre-determined overlap status from the Transcript object
                    overlaps_maincds = getattr(transcript_obj, 'overlaps_maincds', False)
                    logging.debug(f"Using pre-determined overlap status for {transcript_id}: {overlaps_maincds}")

                    # Get codon change
                    codon_change = self.annotator.get_codon_change({
                        'position': variant_coords.transcript,
                        'ref_allele': variant_ref,
                        'alt_allele': variant_alt
                    })
                        
                    # Prepare variant data with all information
                    variant_data = {
                        'position': variant_coords.transcript,
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
                        'transcript_extended': True  # This should always be true for extended transcripts
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
                    display_transcript_id = original_transcript_id
                        
                    result = {
                        'Chromosome': chrom,
                        'Original_Genome_Position': vcf_pos,
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
                        'Transcript_Extended': 'Yes'  # Always mark extended transcripts
                    }
                        
                    results.append(result)
                    logging.debug(f"Successfully processed variant for extended transcript-uORF: {transcript_id}")
            
            # Return all results
            return results
                    
        except Exception as e:
            logging.error(f"Error processing variant: {str(e)}")
            return []
            
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
                logging.debug(f"Using fallback for uORF start: {uorf_start}")
            if uorf_end is None and transcript_obj.uorf_end_genomic is not None:
                uorf_end = transcript_obj.get_transcript_position(transcript_obj.uorf_end_genomic)
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
            'overlaps_maincds': getattr(transcript_obj, 'overlaps_maincds', False)
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