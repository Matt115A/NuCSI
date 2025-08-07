#!/usr/bin/env python3

import os
import re
import gzip
import yaml
import logging
import argparse
import subprocess
import tempfile
from datetime import datetime
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from collections import Counter
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline

# Reverse-complement map
_RC = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(seq):
    return seq.translate(_RC)[::-1]

def load_config(path="configs.yaml"):
    """Load configuration from YAML file."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    
    # Get sequence scanning parameters
    upstream = cfg.get("upstream", "")
    downstream = cfg.get("downstream", "")
    results_base = cfg.get("results_base", "results")
    log_base = os.path.join(results_base, "consensus", "logs")
    
    return upstream, downstream, results_base, log_base

def setup_logging(log_base):
    """Setup logging for the sequence scanning process."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = os.path.join(log_base, timestamp)
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "scan_sequences.log.txt")

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s",
                           datefmt="%Y-%m-%d %H:%M:%S")

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info(f"Logging -> {log_path}")
    return log_dir

def find_sequences_in_read(seq, upstream, downstream, max_length=40):
    """
    Find sequences between upstream and downstream motifs.
    Handles both forward and reverse complement orientations.
    
    Args:
        seq: DNA sequence to search
        upstream: Upstream motif (5'->3')
        downstream: Downstream motif (5'->3')
        max_length: Maximum length of extracted sequence
    
    Returns:
        List of extracted sequences
    """
    sequences = []
    
    # Generate reverse complements of motifs
    upstream_rc = rc(upstream)
    downstream_rc = rc(downstream)
    
    # Search in forward orientation
    sequences.extend(_find_sequences_orientation(seq, upstream, downstream, max_length))
    
    # Search in reverse complement orientation
    sequences.extend(_find_sequences_orientation(seq, downstream_rc, upstream_rc, max_length))
    
    return sequences

def _find_sequences_orientation(seq, upstream, downstream, max_length=40):
    """
    Find sequences between upstream and downstream motifs in a specific orientation.
    Uses fuzzy matching to account for sequencing errors.
    
    Args:
        seq: DNA sequence to search
        upstream: Upstream motif
        downstream: Downstream motif
        max_length: Maximum length of extracted sequence
    
    Returns:
        List of extracted sequences
    """
    sequences = []
    
    # Find all upstream positions with fuzzy matching (allow 1-2 mismatches)
    upstream_positions = []
    for match in re.finditer(upstream, seq, re.IGNORECASE):
        upstream_positions.append(match.end())
    
    # Also try fuzzy matching for upstream
    for i in range(len(seq) - len(upstream) + 1):
        window = seq[i:i+len(upstream)]
        mismatches = sum(1 for a, b in zip(window.upper(), upstream.upper()) if a != b)
        if mismatches <= 2:  # Allow up to 2 mismatches
            upstream_positions.append(i + len(upstream))
    
    # Find all downstream positions with fuzzy matching
    downstream_positions = []
    for match in re.finditer(downstream, seq, re.IGNORECASE):
        downstream_positions.append(match.start())
    
    # Also try fuzzy matching for downstream
    for i in range(len(seq) - len(downstream) + 1):
        window = seq[i:i+len(downstream)]
        mismatches = sum(1 for a, b in zip(window.upper(), downstream.upper()) if a != b)
        if mismatches <= 2:  # Allow up to 2 mismatches
            downstream_positions.append(i)
    
    # Remove duplicates and sort
    upstream_positions = sorted(list(set(upstream_positions)))
    downstream_positions = sorted(list(set(downstream_positions)))
    
    # For each upstream position, find the closest downstream position
    for up_pos in upstream_positions:
        for down_pos in downstream_positions:
            if down_pos > up_pos:  # Downstream must be after upstream
                extracted_seq = seq[up_pos:down_pos]
                if len(extracted_seq) <= max_length and len(extracted_seq) > 0:
                    sequences.append(extracted_seq)
                break  # Take the first downstream position after this upstream
    
    return sequences

def process_fastq_file(file_path, upstream, downstream, max_length=40):
    """
    Process a FASTQ file and extract sequences between motifs.
    
    Args:
        file_path: Path to FASTQ file
        upstream: Upstream motif
        downstream: Downstream motif
        max_length: Maximum sequence length
    
    Returns:
        List of extracted sequences
    """
    sequences = []
    total_reads = 0
    reads_with_upstream = 0
    reads_with_downstream = 0
    reads_with_both = 0
    opener = gzip.open if file_path.endswith(".gz") else open
    
    with opener(file_path, "rt") as fh:
        for record in SeqIO.parse(fh, "fastq"):
            total_reads += 1
            seq = str(record.seq)
            
            # Debug: Count motif occurrences
            upstream_count = len(list(re.finditer(upstream, seq, re.IGNORECASE)))
            downstream_count = len(list(re.finditer(downstream, seq, re.IGNORECASE)))
            upstream_rc_count = len(list(re.finditer(rc(upstream), seq, re.IGNORECASE)))
            downstream_rc_count = len(list(re.finditer(rc(downstream), seq, re.IGNORECASE)))
            
            if upstream_count > 0 or upstream_rc_count > 0:
                reads_with_upstream += 1
            if downstream_count > 0 or downstream_rc_count > 0:
                reads_with_downstream += 1
            if (upstream_count > 0 and downstream_count > 0) or (upstream_rc_count > 0 and downstream_rc_count > 0):
                reads_with_both += 1
            
            extracted = find_sequences_in_read(seq, upstream, downstream, max_length)
            sequences.extend(extracted)
    
    # Log debugging information
    logging.info(f"  Total reads processed: {total_reads}")
    logging.info(f"  Reads with upstream motif: {reads_with_upstream}")
    logging.info(f"  Reads with downstream motif: {reads_with_downstream}")
    logging.info(f"  Reads with both motifs: {reads_with_both}")
    
    return sequences

def run_mafft_alignment(sequences, output_dir):
    """
    Run MAFFT alignment on the sequences.
    
    Args:
        sequences: List of sequences to align
        output_dir: Directory to save results
    
    Returns:
        Path to the aligned FASTA file
    """
    if not sequences:
        logging.warning("No sequences to align")
        return None
    
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
        for i, seq in enumerate(sequences):
            temp_fasta.write(f">sequence_{i:06d}\n{seq}\n")
        temp_fasta_path = temp_fasta.name
    
    # Run MAFFT alignment
    output_fasta = os.path.join(output_dir, "aligned_sequences.fasta")
    mafft_cline = MafftCommandline(input=temp_fasta_path)
    
    try:
        stdout, stderr = mafft_cline()
        with open(output_fasta, 'w') as fh:
            fh.write(stdout)
        logging.info(f"MAFFT alignment completed: {output_fasta}")
    except subprocess.CalledProcessError as e:
        logging.error(f"MAFFT alignment failed: {e}")
        return None
    finally:
        # Clean up temporary file
        os.unlink(temp_fasta_path)
    
    return output_fasta

def create_consensus_logo(aligned_fasta, output_dir):
    """
    Create a consensus logo from the aligned sequences.
    
    Args:
        aligned_fasta: Path to aligned FASTA file
        output_dir: Directory to save results
    """
    if not os.path.exists(aligned_fasta):
        logging.warning("Aligned FASTA file not found")
        return
    
    # Read the alignment
    alignment = AlignIO.read(aligned_fasta, "fasta")
    
    # Convert to position frequency matrix
    seq_length = alignment.get_alignment_length()
    pfm = np.zeros((4, seq_length))  # A, C, G, T
    
    for record in alignment:
        seq = str(record.seq).upper()
        for i, base in enumerate(seq):
            if base == 'A':
                pfm[0, i] += 1
            elif base == 'C':
                pfm[1, i] += 1
            elif base == 'G':
                pfm[2, i] += 1
            elif base == 'T':
                pfm[3, i] += 1
    
    # Convert to probability matrix
    pfm = pfm / pfm.sum(axis=0, keepdims=True)
    
    # Create DataFrame for logomaker
    bases = ['A', 'C', 'G', 'T']
    df = pd.DataFrame(pfm.T, columns=bases)
    
    # Create logo plot
    plt.figure(figsize=(max(8, seq_length/2), 4))
    logomaker.Logo(df, shade_below=0.5, fade_below=0.5)
    plt.title("Consensus Sequence Logo")
    plt.xlabel("Position")
    plt.ylabel("Information (bits)")
    plt.tight_layout()
    
    # Save logo
    logo_path = os.path.join(output_dir, "consensus_logo.png")
    plt.savefig(logo_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Consensus logo saved: {logo_path}")

def main():
    parser = argparse.ArgumentParser(description="Scan sequences between upstream and downstream motifs")
    parser.add_argument("-c", "--config", required=True, help="Path to configs.yaml")
    parser.add_argument("--max-length", type=int, default=40, help="Maximum sequence length (default: 40)")
    args = parser.parse_args()
    
    # Load configuration
    upstream, downstream, results_base, log_base = load_config(args.config)
    log_dir = setup_logging(log_base)
    
    # Create output directory
    output_dir = os.path.join(results_base, "consensus")
    os.makedirs(output_dir, exist_ok=True)
    
    logging.info(f"Upstream motif: {upstream}")
    logging.info(f"Downstream motif: {downstream}")
    logging.info(f"Maximum sequence length: {args.max_length}")
    
    # Find FASTQ files in qc_reads directory
    qc_dir = "inputs/qc_reads"
    fastq_files = []
    for file in os.listdir(qc_dir):
        if file.endswith(('.fastq', '.fastq.gz')) and file.endswith('_qc.fastq.gz'):
            fastq_files.append(os.path.join(qc_dir, file))
    
    if not fastq_files:
        logging.warning("No FASTQ files found in qc_reads directory")
        return
    
    # Process all FASTQ files
    all_sequences = []
    for fastq_file in tqdm(fastq_files, desc="Processing FASTQ files"):
        logging.info(f"Processing {fastq_file}")
        sequences = process_fastq_file(fastq_file, upstream, downstream, args.max_length)
        all_sequences.extend(sequences)
        logging.info(f"Found {len(sequences)} sequences in {os.path.basename(fastq_file)}")
        
        # Debug: Check a few sequences to see what we're finding
        if sequences:
            logging.info(f"Sample sequences from {os.path.basename(fastq_file)}:")
            for i, seq in enumerate(sequences[:3]):
                logging.info(f"  {i+1}: {seq}")
    
    logging.info(f"Total sequences found: {len(all_sequences)}")
    
    if not all_sequences:
        logging.warning("No sequences found matching criteria")
        return
    
    # Save raw sequences
    raw_sequences_path = os.path.join(output_dir, "raw_sequences.fasta")
    with open(raw_sequences_path, 'w') as fh:
        for i, seq in enumerate(all_sequences):
            fh.write(f">sequence_{i:06d}\n{seq}\n")
    
    logging.info(f"Raw sequences saved: {raw_sequences_path}")
    
    # Run MAFFT alignment
    aligned_fasta = run_mafft_alignment(all_sequences, output_dir)
    
    if aligned_fasta:
        # Create consensus logo
        create_consensus_logo(aligned_fasta, output_dir)
        
        # Save summary statistics
        summary_path = os.path.join(output_dir, "summary.txt")
        with open(summary_path, 'w') as fh:
            fh.write(f"Sequence Scanning Summary\n")
            fh.write(f"========================\n")
            fh.write(f"Upstream motif: {upstream}\n")
            fh.write(f"Downstream motif: {downstream}\n")
            fh.write(f"Maximum sequence length: {args.max_length}\n")
            fh.write(f"Total sequences found: {len(all_sequences)}\n")
            fh.write(f"Files processed: {len(fastq_files)}\n")
            fh.write(f"Average sequence length: {np.mean([len(s) for s in all_sequences]):.1f}\n")
            fh.write(f"Sequence length range: {min(len(s) for s in all_sequences)}-{max(len(s) for s in all_sequences)}\n")
        
        logging.info(f"Summary saved: {summary_path}")
    
    logging.info("Sequence scanning completed successfully!")

if __name__ == "__main__":
    main() 