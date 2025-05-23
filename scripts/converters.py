import logging
import gzip
from typing import Dict, Tuple
import pandas as pd

from scripts.models import Exon, Transcript
from scripts.parsers import GTFParser


class CoordinateConverter:
    """Handles conversion between genomic and transcript coordinates."""
    
    def __init__(self, bed_file: str, gtf_file: str, uorf_type: str = "ALL", debug_mode: bool = False):
        """Initialize the converter with input files."""
        self.bed_file = bed_file
        self.gtf_file = gtf_file
        self.transcripts: Dict[str, Transcript] = {}
        self.raw_bed_entries = {}  # Store the raw bed entries for overlap determination
        self.debug_mode = debug_mode
        self.uorf_type = uorf_type  # Store the uORF type option
        
        # Setup logging - use INFO level by default, DEBUG only when debug_mode is True
        log_level = logging.DEBUG if debug_mode else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            force=True
        )
        
        logging.info(f"Initializing coordinate converter with uORF type filter: {uorf_type}")
        
        # Parse files
        self._parse_files()

    def _parse_files(self) -> None:
        """Parse both BED and GTF files to build transcript information."""
        gtf_data = self._parse_gtf()
        self._parse_bed(gtf_data)
        
        # Print summary of processed transcripts
        logging.info(f"Processed {len(self.transcripts)} transcript-uORF pairs")
        logging.info(f"Raw BED entries: {len(self.raw_bed_entries)}")

    def _parse_gtf(self) -> Dict:
        """Parse GTF file to extract CDS coordinates."""
        cds_by_transcript = {}
        exons_by_transcript = {}
        
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parsed = GTFParser.parse_line(line)
                if not parsed:
                    continue

                transcript_id = parsed['transcript_id']
                if not transcript_id:
                    continue

                self._process_gtf_feature(parsed, transcript_id, 
                                       exons_by_transcript, cds_by_transcript)

        return {'cds': cds_by_transcript, 'exons': exons_by_transcript}

    def _process_gtf_feature(self, parsed: dict, transcript_id: str,
                           exons_by_transcript: Dict, cds_by_transcript: Dict) -> None:
        """Process individual features from GTF file."""
        if parsed['feature'] == 'exon':
            if transcript_id not in exons_by_transcript:
                exons_by_transcript[transcript_id] = []
            exon = Exon(
                start=len(exons_by_transcript[transcript_id]),
                length=parsed['end'] - parsed['start'] + 1,
                genome_start=parsed['start'],
                genome_end=parsed['end']
            )
            exons_by_transcript[transcript_id].append(exon)
        
        elif parsed['feature'] == 'CDS':
            if transcript_id not in cds_by_transcript:
                cds_by_transcript[transcript_id] = {'starts': [], 'ends': []}
            cds_by_transcript[transcript_id]['starts'].append(parsed['start'])
            cds_by_transcript[transcript_id]['ends'].append(parsed['end'])

    def _parse_bed(self, gtf_data: Dict) -> None:
        """Parse BED file and create transcript objects, applying uORF type filter if specified."""
        bed_df = pd.read_csv(self.bed_file, sep='\t', header=None)
        
        logging.info(f"BED file contains {len(bed_df)} entries")
        
        # Create a dictionary to track which transcripts have been processed
        processed_transcripts = {}
        skipped_frames = 0
        processed_frames = 0
        
        for _, bed_row in bed_df.iterrows():
            transcript_id = bed_row[3].split('|')[0].split('.')[0]
            
            # Store raw bed entry before any processing
            entry = {
                'chromosome': bed_row[0],
                'start': int(bed_row[1]),
                'end': int(bed_row[2]),
                'name': bed_row[3],
                'strand': bed_row[5],
                'transcript_id': transcript_id
            }
            
            # Check if this entry has start codon information in the name field
            # BED name format could include frame type: "transcript_id|other_info|ATG" or "transcript_id|other_info|non-ATG"
            has_atg_start = False
            is_non_atg = False
            
            # Try to extract frame type from the name field if available
            bed_name_parts = bed_row[3].split('|')
            if len(bed_name_parts) >= 3:
                # Check the last part for ATG indicator
                frame_indicator = bed_name_parts[-1].upper()
                if frame_indicator == "ATG":
                    has_atg_start = True
                elif frame_indicator in ["NON-ATG", "NONATG", "NON_ATG"]:
                    is_non_atg = True
            
            # Apply uORF type filter if specified
            if self.uorf_type != "ALL":
                if self.uorf_type == "ATG" and not has_atg_start:
                    skipped_frames += 1
                    continue
                elif self.uorf_type == "NON-ATG" and not is_non_atg:
                    skipped_frames += 1
                    continue
            
            processed_frames += 1
            
            # Create a unique key for this uORF
            uorf_id = f"{transcript_id}_uorf_{len(processed_transcripts.get(transcript_id, []))+1}"
            
            # Store raw bed entry with its own ID to preserve uniqueness
            self.raw_bed_entries[uorf_id] = entry
            
            if transcript_id in gtf_data['exons']:
                mainorf_coords = self._get_mainorf_coords(transcript_id, gtf_data['cds'])
                
                # Log what we're processing
                logging.info(f"Processing BED row: {entry}")
                logging.info(f"Creating transcript {uorf_id} with uORF at {entry['start']+1}-{entry['end']}")

                # Always create a new transcript with unique ID for each uORF
                new_transcript = Transcript(
                    transcript_id=uorf_id,
                    chromosome=entry['chromosome'],
                    strand=entry['strand'],
                    exons=gtf_data['exons'][transcript_id].copy() if transcript_id in gtf_data['exons'] else [],
                    mainorf_start=mainorf_coords[0],
                    mainorf_end=mainorf_coords[1],
                    uorf_start=int(bed_row[1]) + 1,  # BED is 0-based, convert to 1-based
                    uorf_end=int(bed_row[2]),
                    start_codon_type="ATG" if has_atg_start else "NON-ATG" if is_non_atg else None,  # Store the start codon type
                    debug_mode=self.debug_mode  # Pass debug_mode to Transcript
                )
                
                # Store the new transcript
                self.transcripts[uorf_id] = new_transcript
                
                # Track that we've processed this transcript
                if transcript_id not in processed_transcripts:
                    processed_transcripts[transcript_id] = []
                processed_transcripts[transcript_id].append(uorf_id)
            else:
                logging.warning(f"No exons found for transcript {transcript_id}")

        # Log summary of uORF type filtering
        if self.uorf_type != "ALL":
            logging.info(f"uORF type filter '{self.uorf_type}' applied: processed {processed_frames} entries, skipped {skipped_frames} entries")
            
        # Log summary of transcripts with multiple uORFs
        for transcript_id, uorf_ids in processed_transcripts.items():
            if len(uorf_ids) > 1:
                logging.info(f"Transcript {transcript_id} has {len(uorf_ids)} uORFs: {', '.join(uorf_ids)}")
                
        # After all transcripts are created, determine overlaps using raw coordinates
        self._determine_overlaps()

    def _determine_overlaps(self):
        """
        Determine which uORFs overlap with mainCDS based on raw genomic coordinates
        from the BED and GTF files.
        """
        logging.info("Determining uORF-mainCDS overlaps...")
        for transcript_id, transcript in self.transcripts.items():
            # Get the raw entry
            raw_entry = self.raw_bed_entries.get(transcript_id)
            if not raw_entry:
                logging.warning(f"No raw BED entry found for {transcript_id}")
                continue
                
            # A uORF and mainCDS overlap if:
            # For + strand: uORF end >= mainCDS start
            # For - strand: uORF start <= mainCDS end
            if transcript.strand == '+':
                # Raw uORF end is in BED format (0-based, exclusive)
                uorf_end = raw_entry['end']
                # MainORF start is in GTF format (1-based, inclusive)
                mainorf_start = transcript.mainorf_start_genomic
                
                if mainorf_start is not None:
                    transcript.overlaps_maincds = (uorf_end >= mainorf_start)
                    logging.info(f"Transcript {transcript_id}: uORF end={uorf_end}, mainCDS start={mainorf_start}, overlap={transcript.overlaps_maincds}")
            else:
                # For - strand, coordinates might be swapped in the raw data
                uorf_start = raw_entry['start'] + 1  # Adjust to 1-based
                mainorf_end = transcript.mainorf_end_genomic
                
                if mainorf_end is not None:
                    transcript.overlaps_maincds = (uorf_start <= mainorf_end)
                    logging.info(f"Transcript {transcript_id}: uORF start={uorf_start}, mainCDS end={mainorf_end}, overlap={transcript.overlaps_maincds}")

    def _get_mainorf_coords(self, transcript_id: str, cds_data: Dict) -> Tuple[int, int]:
        """Get main ORF coordinates for a transcript."""
        if transcript_id in cds_data:
            return (
                min(cds_data[transcript_id]['starts']),
                max(cds_data[transcript_id]['ends'])
            )
        return None, None

    def genome_to_transcript_pos(self, transcript_id: str, genome_pos: int) -> str:
        """Convert genomic position to transcript position."""
        transcript = self.transcripts.get(transcript_id)
        if not transcript:
            logging.warning(f"Transcript {transcript_id} not found")
            return "NA"

        return transcript.genome_to_transcript.get(genome_pos, "NA")