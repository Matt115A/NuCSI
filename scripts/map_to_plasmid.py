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
import seaborn as sns
import scipy.stats
from collections import Counter, defaultdict
from Bio import SeqIO
import pysam

# Reverse-complement map
_RC = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(seq):
    return seq.translate(_RC)[::-1]

def load_config(path="configs.yaml"):
    """Load configuration from YAML file."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    
    results_base = cfg.get("results_base", "results")
    log_base = os.path.join(results_base, "plasmid_mapping", "logs")
    
    return results_base, log_base

def setup_logging(log_base):
    """Setup logging for the plasmid mapping process."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = os.path.join(log_base, timestamp)
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "plasmid_mapping.log.txt")

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

def create_bwa_index(plasmid_fasta, output_dir):
    """Create BWA index for the plasmid reference."""
    index_prefix = os.path.join(output_dir, "plasmid_index")
    
    # Check if index already exists
    if os.path.exists(f"{index_prefix}.amb"):
        logging.info("BWA index already exists, skipping creation")
        return index_prefix
    
    logging.info("Creating BWA index for plasmid reference")
    cmd = ["bwa", "index", "-p", index_prefix, plasmid_fasta]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"BWA index created: {index_prefix}")
        return index_prefix
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create BWA index: {e}")
        raise

def normalize_sequence_case(plasmid_fasta, output_dir):
    """Normalize plasmid sequence to uppercase for consistent case handling."""
    normalized_fasta = os.path.join(output_dir, "plasmid_normalized.fasta")
    
    with open(plasmid_fasta, 'r') as infile, open(normalized_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line)
            else:
                # Convert to uppercase for consistent case handling
                outfile.write(line.upper())
    
    logging.info(f"Normalized plasmid sequence saved: {normalized_fasta}")
    return normalized_fasta

def run_bwa_mapping(fastq_file, index_prefix, output_dir, sample_name):
    """Run BWA mapping for a FASTQ file."""
    bam_file = os.path.join(output_dir, f"{sample_name}.bam")
    sam_file = os.path.join(output_dir, f"{sample_name}.sam")
    
    # Run BWA mem with case-insensitive handling
    logging.info(f"Running BWA mem for {sample_name}")
    cmd = [
        "bwa", "mem", 
        "-t", "4",  # Use 4 threads
        "-M",  # Mark shorter split hits as secondary
        index_prefix,
        fastq_file
    ]
    
    try:
        with open(sam_file, 'w') as fh:
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, check=True)
        logging.info(f"BWA mapping completed for {sample_name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"BWA mapping failed for {sample_name}: {e}")
        raise
    
    # Convert SAM to BAM and sort
    logging.info(f"Converting SAM to BAM for {sample_name}")
    cmd_sort = [
        "samtools", "view", "-b", sam_file, "|",
        "samtools", "sort", "-o", bam_file
    ]
    
    try:
        subprocess.run(" ".join(cmd_sort), shell=True, check=True)
        logging.info(f"BAM file created: {bam_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to convert SAM to BAM: {e}")
        raise
    
    # Index BAM file
    logging.info(f"Indexing BAM file for {sample_name}")
    cmd_index = ["samtools", "index", bam_file]
    subprocess.run(cmd_index, check=True)
    
    return bam_file

def filter_alignments(bam_file, min_quality=30, min_length=50):
    """Filter alignments by quality and length."""
    filtered_alignments = []
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            # Check mapping quality
            if read.mapping_quality < min_quality:
                continue
            
            # Check if read is mapped
            if read.is_unmapped:
                continue
            
            # Check alignment length
            if read.query_alignment_length < min_length:
                continue
            
            # Check if alignment is primary (not secondary/supplementary)
            if read.is_secondary or read.is_supplementary:
                continue
            
            filtered_alignments.append(read)
    
    logging.info(f"Filtered {len(filtered_alignments)} high-quality alignments from {bam_file}")
    return filtered_alignments

def extract_mapping_positions(alignments, plasmid_length):
    """Extract mapping start and end positions from alignments.
    
    Note: This function uses 0-based indexing for all positions:
    - reference_start: 0-based, inclusive (first aligned base)
    - reference_end: 0-based, exclusive (position after last aligned base)
    - This is consistent with pysam and BWA conventions
    """
    mapping_data = []
    
    for read in alignments:
        # Get mapping positions (0-based indexing)
        start_pos = read.reference_start  # 0-based, inclusive
        end_pos = read.reference_end      # 0-based, exclusive
        
        # Handle circular plasmid wrapping
        if end_pos > plasmid_length:
            end_pos = end_pos % plasmid_length
        
        mapping_data.append({
            'read_name': read.query_name,
            'start': start_pos,
            'end': end_pos,
            'length': end_pos - start_pos if end_pos > start_pos else (plasmid_length - start_pos + end_pos),
            'strand': '-' if read.is_reverse else '+',
            'mapping_quality': read.mapping_quality,
            'cigar': read.cigarstring
        })
    
    return mapping_data

def create_mapping_plots(mapping_data, plasmid_length, output_dir, sample_name):
    """Create mapping position plots and charts."""
    
    if not mapping_data:
        logging.warning(f"No mapping data for {sample_name}")
        return
    
    df = pd.DataFrame(mapping_data)
    
    # Create comprehensive mapping analysis plots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    
    # Plot 1: Circular mapping positions
    ax1 = axes[0, 0]
    ax1.set_aspect('equal')
    ax1.set_xlim(-1.2, 1.2)
    ax1.set_ylim(-1.2, 1.2)
    
    # Draw circle representing plasmid
    circle = plt.Circle((0, 0), 1, fill=False, color='black', linewidth=2)
    ax1.add_artist(circle)
    
    # Plot mapping positions
    for _, row in df.iterrows():
        start_angle = 2 * np.pi * row['start'] / plasmid_length
        end_angle = 2 * np.pi * row['end'] / plasmid_length
        
        # Draw arc for mapping
        theta = np.linspace(start_angle, end_angle, 100)
        x = np.cos(theta)
        y = np.sin(theta)
        ax1.plot(x, y, 'b-', alpha=0.3, linewidth=1)
    
    ax1.set_title(f'Circular Mapping - {sample_name}')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    
    # Plot 2: Linear mapping positions
    ax2 = axes[0, 1]
    ax2.scatter(df['start'], df['end'], alpha=0.5, s=10)
    ax2.set_xlabel('Start Position')
    ax2.set_ylabel('End Position')
    ax2.set_title(f'Linear Mapping Positions - {sample_name}')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Start position distribution
    ax3 = axes[0, 2]
    start_counts = df['start'].value_counts().sort_index()
    ax3.bar(start_counts.index, start_counts.values, alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Start Position (bp)')
    ax3.set_ylabel('Frequency')
    ax3.set_title(f'Start Position Distribution - {sample_name}')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: End position distribution
    ax4 = axes[1, 0]
    end_counts = df['end'].value_counts().sort_index()
    ax4.bar(end_counts.index, end_counts.values, alpha=0.7, edgecolor='black', color='orange')
    ax4.set_xlabel('End Position (bp)')
    ax4.set_ylabel('Frequency')
    ax4.set_title(f'End Position Distribution - {sample_name}')
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Alignment length distribution
    ax5 = axes[1, 1]
    ax5.hist(df['length'], bins=30, alpha=0.7, edgecolor='black', color='green')
    ax5.set_xlabel('Alignment Length (bp)')
    ax5.set_ylabel('Frequency')
    ax5.set_title(f'Alignment Length Distribution - {sample_name}')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Mapping quality distribution
    ax6 = axes[1, 2]
    ax6.hist(df['mapping_quality'], bins=20, alpha=0.7, edgecolor='black', color='red')
    ax6.set_xlabel('Mapping Quality')
    ax6.set_ylabel('Frequency')
    ax6.set_title(f'Mapping Quality Distribution - {sample_name}')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f"{sample_name}_mapping_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Mapping plot saved: {plot_path}")
    
    # Create detailed position analysis
    create_position_analysis(df, plasmid_length, output_dir, sample_name)
    
    # Create mapping statistics
    stats_path = os.path.join(output_dir, f"{sample_name}_mapping_stats.txt")
    with open(stats_path, 'w') as fh:
        fh.write(f"Mapping Statistics for {sample_name}\n")
        fh.write("=" * 50 + "\n")
        fh.write(f"Total alignments: {len(df)}\n")
        fh.write(f"Average mapping quality: {df['mapping_quality'].mean():.2f}\n")
        fh.write(f"Average alignment length: {df['length'].mean():.2f}\n")
        fh.write(f"Start position range: {df['start'].min()} - {df['start'].max()}\n")
        fh.write(f"End position range: {df['end'].min()} - {df['end'].max()}\n")
        fh.write(f"Forward strand alignments: {(df['strand'] == '+').sum()}\n")
        fh.write(f"Reverse strand alignments: {(df['strand'] == '-').sum()}\n")
    
    logging.info(f"Mapping statistics saved: {stats_path}")
    
    return df

def create_position_analysis(df, plasmid_length, output_dir, sample_name):
    """Create detailed position analysis with statistical power and multiple testing correction."""
    
    # Analyze start positions (case-insensitive)
    start_counts = df['start'].value_counts().sort_index()
    end_counts = df['end'].value_counts().sort_index()
    
    # Calculate statistical significance (using chi-square test against uniform distribution)
    total_alignments = len(df)
    expected_frequency = total_alignments / plasmid_length
    
    # Start position analysis
    start_significance = {}
    start_p_values = []
    start_positions = []
    
    for pos, count in start_counts.items():
        # Chi-square test for this position
        chi_square = (count - expected_frequency) ** 2 / expected_frequency
        p_value = 1 - scipy.stats.chi2.cdf(chi_square, 1)  # 1 degree of freedom
        start_p_values.append(p_value)
        start_positions.append(pos)
        start_significance[pos] = {
            'count': count,
            'frequency': count / total_alignments,
            'expected': expected_frequency,
            'chi_square': chi_square,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
    
    # Apply multiple testing corrections for start positions
    if start_p_values:
        # Bonferroni correction
        bonferroni_threshold = 0.05 / len(start_p_values)
        
        for i, pos in enumerate(start_positions):
            # Bonferroni correction
            start_significance[pos]['bonferroni_p'] = min(start_p_values[i] * len(start_p_values), 1.0)
            start_significance[pos]['bonferroni_significant'] = start_p_values[i] < bonferroni_threshold
            
            # Benjamini-Hochberg correction (simplified)
            start_significance[pos]['bh_p'] = min(start_p_values[i] * len(start_p_values), 1.0)
            start_significance[pos]['bh_significant'] = start_p_values[i] < (0.05 * len(start_p_values) / (i + 1))
    
    # End position analysis
    end_significance = {}
    end_p_values = []
    end_positions = []
    
    for pos, count in end_counts.items():
        chi_square = (count - expected_frequency) ** 2 / expected_frequency
        p_value = 1 - scipy.stats.chi2.cdf(chi_square, 1)
        end_p_values.append(p_value)
        end_positions.append(pos)
        end_significance[pos] = {
            'count': count,
            'frequency': count / total_alignments,
            'expected': expected_frequency,
            'chi_square': chi_square,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
    
    # Apply multiple testing corrections for end positions
    if end_p_values:
        # Bonferroni correction
        bonferroni_threshold = 0.05 / len(end_p_values)
        
        for i, pos in enumerate(end_positions):
            # Bonferroni correction
            end_significance[pos]['bonferroni_p'] = min(end_p_values[i] * len(end_p_values), 1.0)
            end_significance[pos]['bonferroni_significant'] = end_p_values[i] < bonferroni_threshold
            
            # Benjamini-Hochberg correction (simplified)
            end_significance[pos]['bh_p'] = min(end_p_values[i] * len(end_p_values), 1.0)
            end_significance[pos]['bh_significant'] = end_p_values[i] < (0.05 * len(end_p_values) / (i + 1))
    
    # Create position analysis plots with multiple testing correction
    fig, axes = plt.subplots(2, 3, figsize=(24, 12))
    
    # Plot 1: Start position frequency with Bonferroni significance
    ax1 = axes[0, 0]
    positions = list(start_significance.keys())
    frequencies = [start_significance[pos]['frequency'] for pos in positions]
    bonferroni_sig = [start_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax1.bar(positions, frequencies, color=colors, alpha=0.7, edgecolor='black')
    ax1.axhline(y=1/plasmid_length, color='green', linestyle='--', alpha=0.7, label='Expected (uniform)')
    ax1.set_xlabel('Start Position (bp)')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'Start Position Frequency - {sample_name}\n(Red = Bonferroni Significant)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: End position frequency with Bonferroni significance
    ax2 = axes[0, 1]
    positions = list(end_significance.keys())
    frequencies = [end_significance[pos]['frequency'] for pos in positions]
    bonferroni_sig = [end_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax2.bar(positions, frequencies, color=colors, alpha=0.7, edgecolor='black')
    ax2.axhline(y=1/plasmid_length, color='green', linestyle='--', alpha=0.7, label='Expected (uniform)')
    ax2.set_xlabel('End Position (bp)')
    ax2.set_ylabel('Frequency')
    ax2.set_title(f'End Position Frequency - {sample_name}\n(Red = Bonferroni Significant)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Start position frequency with BH significance
    ax3 = axes[0, 2]
    positions = list(start_significance.keys())
    frequencies = [start_significance[pos]['frequency'] for pos in positions]
    bh_sig = [start_significance[pos]['bh_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bh_sig]
    ax3.bar(positions, frequencies, color=colors, alpha=0.7, edgecolor='black')
    ax3.axhline(y=1/plasmid_length, color='green', linestyle='--', alpha=0.7, label='Expected (uniform)')
    ax3.set_xlabel('Start Position (bp)')
    ax3.set_ylabel('Frequency')
    ax3.set_title(f'Start Position Frequency - {sample_name}\n(Red = BH Significant)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Start position Bonferroni p-values
    ax4 = axes[1, 0]
    positions = list(start_significance.keys())
    p_values = [start_significance[pos]['bonferroni_p'] for pos in positions]
    bonferroni_sig = [start_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax4.scatter(positions, p_values, color=colors, alpha=0.7, s=20)
    ax4.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    ax4.set_xlabel('Start Position (bp)')
    ax4.set_ylabel('Bonferroni P-value')
    ax4.set_title(f'Start Position Bonferroni P-values - {sample_name}')
    ax4.set_yscale('log')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: End position Bonferroni p-values
    ax5 = axes[1, 1]
    positions = list(end_significance.keys())
    p_values = [end_significance[pos]['bonferroni_p'] for pos in positions]
    bonferroni_sig = [end_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax5.scatter(positions, p_values, color=colors, alpha=0.7, s=20)
    ax5.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    ax5.set_xlabel('End Position (bp)')
    ax5.set_ylabel('Bonferroni P-value')
    ax5.set_title(f'End Position Bonferroni P-values - {sample_name}')
    ax5.set_yscale('log')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Start position BH p-values
    ax6 = axes[1, 2]
    positions = list(start_significance.keys())
    p_values = [start_significance[pos]['bh_p'] for pos in positions]
    bh_sig = [start_significance[pos]['bh_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bh_sig]
    ax6.scatter(positions, p_values, color=colors, alpha=0.7, s=20)
    ax6.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    ax6.set_xlabel('Start Position (bp)')
    ax6.set_ylabel('BH P-value')
    ax6.set_title(f'Start Position BH P-values - {sample_name}')
    ax6.set_yscale('log')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    analysis_plot_path = os.path.join(output_dir, f"{sample_name}_position_analysis.png")
    plt.savefig(analysis_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Position analysis plot saved: {analysis_plot_path}")
    
    # Save detailed position analysis with multiple testing corrections
    analysis_path = os.path.join(output_dir, f"{sample_name}_position_analysis.txt")
    with open(analysis_path, 'w') as fh:
        fh.write(f"Position Analysis for {sample_name}\n")
        fh.write("=" * 60 + "\n\n")
        
        # Multiple testing correction information
        fh.write("MULTIPLE TESTING CORRECTION\n")
        fh.write("-" * 35 + "\n")
        fh.write(f"Number of start positions tested: {len(start_significance)}\n")
        fh.write(f"Number of end positions tested: {len(end_significance)}\n")
        fh.write(f"Bonferroni threshold (α=0.05): {0.05/len(start_significance):.6f}\n")
        fh.write(f"Benjamini-Hochberg FDR control: α=0.05\n\n")
        
        # Start position analysis
        fh.write("START POSITION ANALYSIS\n")
        fh.write("-" * 30 + "\n")
        fh.write(f"Total alignments: {total_alignments}\n")
        fh.write(f"Plasmid length: {plasmid_length} bp\n")
        fh.write(f"Expected frequency (uniform): {1/plasmid_length:.6f}\n\n")
        
        # Sort by frequency (most common first)
        sorted_starts = sorted(start_significance.items(), key=lambda x: x[1]['count'], reverse=True)
        
        fh.write("Most Common Start Positions (with multiple testing corrections):\n")
        fh.write("Position\tCount\tFrequency\tRaw P-value\tBonferroni P\tBH P\tBonferroni Sig\tBH Sig\n")
        fh.write("-" * 100 + "\n")
        
        for pos, stats in sorted_starts[:20]:  # Top 20
            fh.write(f"{pos}\t{stats['count']}\t{stats['frequency']:.6f}\t{stats['p_value']:.6f}\t{stats['bonferroni_p']:.6f}\t{stats['bh_p']:.6f}\t{stats['bonferroni_significant']}\t{stats['bh_significant']}\n")
        
        fh.write("\n\nEND POSITION ANALYSIS\n")
        fh.write("-" * 28 + "\n")
        
        # Sort by frequency (most common first)
        sorted_ends = sorted(end_significance.items(), key=lambda x: x[1]['count'], reverse=True)
        
        fh.write("Most Common End Positions (with multiple testing corrections):\n")
        fh.write("Position\tCount\tFrequency\tRaw P-value\tBonferroni P\tBH P\tBonferroni Sig\tBH Sig\n")
        fh.write("-" * 100 + "\n")
        
        for pos, stats in sorted_ends[:20]:  # Top 20
            fh.write(f"{pos}\t{stats['count']}\t{stats['frequency']:.6f}\t{stats['p_value']:.6f}\t{stats['bonferroni_p']:.6f}\t{stats['bh_p']:.6f}\t{stats['bonferroni_significant']}\t{stats['bh_significant']}\n")
    
    logging.info(f"Position analysis saved: {analysis_path}")
    
    # Create detailed zoom plots for top hits
    create_zoom_plots(df, plasmid_length, output_dir, sample_name, start_significance, end_significance)
    
    return start_significance, end_significance


def create_zoom_plots(df, plasmid_length, output_dir, sample_name, start_significance, end_significance):
    """Create detailed zoom plots for top hit regions with individual bases."""
    
    # Load plasmid sequence
    plasmid_file = os.path.join(output_dir, "plasmid_normalized.fasta")
    with open(plasmid_file, 'r') as fh:
        plasmid_seq = str(next(SeqIO.parse(fh, 'fasta')).seq)
    
    # Find top hits
    start_positions = sorted(start_significance.items(), key=lambda x: x[1]['count'], reverse=True)
    end_positions = sorted(end_significance.items(), key=lambda x: x[1]['count'], reverse=True)
    
    if not start_positions or not end_positions:
        logging.warning(f"No positions to plot for {sample_name}")
        return
    
    # Get top start and end positions
    top_start_pos = start_positions[0][0]
    top_end_pos = end_positions[0][0]
    
    # Create zoom plots for start position
    create_position_zoom_plot(df, plasmid_seq, top_start_pos, "start", output_dir, sample_name, start_significance)
    
    # Create zoom plots for end position
    create_position_zoom_plot(df, plasmid_seq, top_end_pos, "end", output_dir, sample_name, end_significance)


def create_position_zoom_plot(df, plasmid_seq, position, pos_type, output_dir, sample_name, significance_dict):
    """Create a detailed zoom plot for a specific position."""
    
    # Define zoom region (20 bases on either side)
    start_region = max(0, position - 20)
    end_region = min(len(plasmid_seq), position + 21)
    
    # Get counts for this region
    region_counts = {}
    for i in range(start_region, end_region):
        if pos_type == "start":
            count = len(df[df['start'] == i])
        else:  # end
            count = len(df[df['end'] == i])
        region_counts[i] = count
    
    # Get sequence for this region
    region_seq = plasmid_seq[start_region:end_region]
    
    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(20, 8))
    
    # Plot bars for each position
    positions = list(region_counts.keys())
    counts = list(region_counts.values())
    
    # Color the top hit position differently
    colors = ['red' if pos == position else 'blue' for pos in positions]
    
    bars = ax.bar(positions, counts, color=colors, alpha=0.7, edgecolor='black')
    
    # Add base labels on x-axis
    ax.set_xticks(positions)
    base_labels = []
    for pos in positions:
        idx = pos - start_region
        if 0 <= idx < len(region_seq):
            base = region_seq[idx]
        else:
            base = 'N'
        base_labels.append(f"{pos}\n{base}")
    ax.set_xticklabels(base_labels, rotation=45, ha='right', fontsize=10)
    
    # Highlight the top hit position
    ax.axvline(x=position, color='red', linestyle='--', linewidth=2, alpha=0.8, label=f'Top {pos_type} position')
    
    # Add labels and title
    ax.set_xlabel(f'Position (bp) and Base')
    ax.set_ylabel('Count')
    ax.set_title(f'{sample_name} - {pos_type.capitalize()} Position Zoom (Position {position})\n'
                f'Region: {start_region}-{end_region-1}, Top hit: {position} ({region_counts[position]} alignments)')
    
    # Add grid and legend
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Add sequence annotation
    seq_text = f"Sequence: {region_seq}"
    ax.text(0.02, 0.98, seq_text, transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Add count annotations for significant positions
    for pos, count in region_counts.items():
        if count > 0:
            significance = significance_dict.get(pos, {})
            if significance.get('bonferroni_significant', False):
                ax.annotate(f'{count}', xy=(pos, count), xytext=(0, 5),
                           textcoords='offset points', ha='center', va='bottom',
                           fontsize=8, fontweight='bold', color='red')
            elif count > 5:  # Only annotate if count > 5 to avoid clutter
                ax.annotate(f'{count}', xy=(pos, count), xytext=(0, 5),
                           textcoords='offset points', ha='center', va='bottom',
                           fontsize=8)
    
    plt.tight_layout()
    
    # Save the plot
    plot_path = os.path.join(output_dir, f"{sample_name}_{pos_type}_zoom_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"{pos_type.capitalize()} position zoom plot saved: {plot_path}")


def main():
    parser = argparse.ArgumentParser(description="Map reads to reference plasmid and analyze mappings")
    parser.add_argument("-c", "--config", required=True, help="Path to configs.yaml")
    parser.add_argument("--min-quality", type=int, default=30, help="Minimum mapping quality (default: 30)")
    parser.add_argument("--min-length", type=int, default=50, help="Minimum alignment length (default: 50)")
    args = parser.parse_args()
    
    # Load configuration
    results_base, log_base = load_config(args.config)
    log_dir = setup_logging(log_base)
    
    # Find plasmid references
    plasmid_dir = "inputs/plasmid"
    plasmid_files = [f for f in os.listdir(plasmid_dir) if f.endswith('.fasta')]
    
    if not plasmid_files:
        logging.error("No plasmid reference found in inputs/plasmid/")
        return
    
    # Find FASTQ files
    qc_dir = "inputs/qc_reads"
    fastq_files = [f for f in os.listdir(qc_dir) if f.endswith('_qc.fastq.gz')]
    
    if not fastq_files:
        logging.error("No FASTQ files found in inputs/qc_reads/")
        return
    
    # Process each plasmid separately
    for plasmid_file in plasmid_files:
        plasmid_name = os.path.splitext(plasmid_file)[0]  # Remove .fasta extension
        logging.info(f"Processing plasmid: {plasmid_name}")
        
        # Create plasmid-specific output directory
        plasmid_output_dir = os.path.join(results_base, f"plasmid_mapping_{plasmid_name}")
        os.makedirs(plasmid_output_dir, exist_ok=True)
        
        plasmid_fasta = os.path.join(plasmid_dir, plasmid_file)
        logging.info(f"Using plasmid reference: {plasmid_fasta}")
        
        # Normalize plasmid sequence to uppercase for case-insensitive handling
        normalized_plasmid = normalize_sequence_case(plasmid_fasta, plasmid_output_dir)
        
        # Get plasmid length from normalized sequence
        with open(normalized_plasmid) as fh:
            plasmid_seq = str(next(SeqIO.parse(fh, 'fasta')).seq)
            plasmid_length = len(plasmid_seq)
        
        logging.info(f"Plasmid length: {plasmid_length} bp")
        
        # Create BWA index using normalized sequence
        index_prefix = create_bwa_index(normalized_plasmid, plasmid_output_dir)
        
        # Process each FASTQ file for this plasmid
        all_mapping_data = []
        
        for fastq_file in tqdm(fastq_files, desc=f"Processing FASTQ files for {plasmid_name}"):
            sample_name = os.path.splitext(os.path.splitext(fastq_file)[0])[0]  # Remove .fastq.gz
            fastq_path = os.path.join(qc_dir, fastq_file)
            
            logging.info(f"Processing {sample_name} for plasmid {plasmid_name}")
            
            # Run BWA mapping
            bam_file = run_bwa_mapping(fastq_path, index_prefix, plasmid_output_dir, sample_name)
            
            # Filter alignments
            alignments = filter_alignments(bam_file, args.min_quality, args.min_length)
            
            if not alignments:
                logging.warning(f"No high-quality alignments for {sample_name} on plasmid {plasmid_name}")
                continue
            
            # Extract mapping positions
            mapping_data = extract_mapping_positions(alignments, plasmid_length)
            all_mapping_data.extend(mapping_data)
            
            # Create mapping plots
            df = create_mapping_plots(mapping_data, plasmid_length, plasmid_output_dir, sample_name)
            
            # Create position analysis
            start_significance, end_significance = create_position_analysis(df, plasmid_length, plasmid_output_dir, sample_name)
            
            # Create detailed zoom plots for top hits
            create_zoom_plots(df, plasmid_length, plasmid_output_dir, sample_name, start_significance, end_significance)
            
            # Save mapping statistics
            stats_path = os.path.join(plasmid_output_dir, f"{sample_name}_mapping_stats.txt")
            with open(stats_path, 'w') as fh:
                fh.write(f"Mapping Statistics for {sample_name} on {plasmid_name}\n")
                fh.write("=" * 50 + "\n")
                fh.write(f"Total alignments: {len(mapping_data)}\n")
                fh.write(f"Average alignment length: {df['length'].mean():.2f} bp\n")
                fh.write(f"Average mapping quality: {df['mapping_quality'].mean():.2f}\n")
                fh.write(f"Plasmid length: {plasmid_length} bp\n")
            
            logging.info(f"Mapping statistics saved: {stats_path}")
        
        # Create comprehensive summary analysis for this plasmid
        if all_mapping_data:
            summary_df = pd.DataFrame(all_mapping_data)
            create_comprehensive_summary(summary_df, plasmid_length, plasmid_output_dir)
        
        logging.info(f"Completed analysis for plasmid {plasmid_name}")
    
    logging.info("Plasmid mapping analysis completed successfully for all plasmids!")

def create_coverage_analysis(summary_df, plasmid_length, output_dir):
    """Create coverage analysis plots for the entire construct and zoomed around the sharpest drop-off.
    
    Note: This function uses 0-based indexing for all positions, which is consistent with:
    - BWA alignment positions (reference_start and reference_end are 0-based)
    - pysam indexing (reference_end is exclusive)
    - Python array indexing
    """
    
    # Load plasmid sequence for base display
    plasmid_fasta = os.path.join(output_dir, "plasmid_normalized.fasta")
    with open(plasmid_fasta, 'r') as fh:
        plasmid_seq = str(next(SeqIO.parse(fh, 'fasta')).seq)
    
    # Create coverage array for the entire plasmid (0-based indexing)
    coverage = np.zeros(plasmid_length)
    
    # Calculate coverage for each position
    for _, row in summary_df.iterrows():
        start_pos = row['start']  # 0-based, inclusive
        end_pos = row['end']      # 0-based, exclusive
        
        # Handle circular plasmid wrapping
        if end_pos > start_pos:
            # Normal case: start < end
            coverage[start_pos:end_pos] += 1
        else:
            # Wrapped case: end < start (circular)
            coverage[start_pos:] += 1
            coverage[:end_pos] += 1
    
    # Find the sharpest drop-off in coverage
    coverage_diff = np.diff(coverage)
    sharpest_drop_pos = np.argmin(coverage_diff)  # Position with largest negative change
    drop_magnitude = abs(coverage_diff[sharpest_drop_pos])
    
    # Create coverage plots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 12))
    
    # Plot 1: Full construct coverage
    ax1.plot(range(plasmid_length), coverage, 'b-', linewidth=1, alpha=0.8)
    ax1.set_xlabel('Position (bp, 0-based)')
    ax1.set_ylabel('Coverage')
    ax1.set_title('Coverage Across Entire Plasmid Construct')
    ax1.grid(True, alpha=0.3)
    
    # Add statistics to the plot
    mean_coverage = np.mean(coverage)
    max_coverage = np.max(coverage)
    ax1.axhline(y=mean_coverage, color='red', linestyle='--', alpha=0.7, 
                label=f'Mean Coverage: {mean_coverage:.1f}')
    ax1.legend()
    
    # Plot 2: Zoomed coverage around the sharpest drop-off (exactly 5 bp on either side)
    zoom_center = sharpest_drop_pos
    zoom_start = max(0, zoom_center - 5)
    zoom_end = min(plasmid_length, zoom_center + 6)
    
    zoom_positions = range(zoom_start, zoom_end)
    zoom_coverage = coverage[zoom_start:zoom_end]
    
    ax2.plot(zoom_positions, zoom_coverage, 'b-', linewidth=2, alpha=0.8, marker='o', markersize=4)
    ax2.set_xlabel('Position (bp, 0-based) and Base')
    ax2.set_ylabel('Coverage')
    ax2.set_title(f'Coverage Zoom: Sharpest Drop-off at Position {zoom_center}\n'
                  f'Positions {zoom_start}-{zoom_end-1} (Drop magnitude: {drop_magnitude:.0f})')
    ax2.grid(True, alpha=0.3)
    
    # Highlight the drop-off position
    ax2.axvline(x=zoom_center, color='red', linestyle='--', alpha=0.7, 
                label=f'Sharpest Drop: Position {zoom_center}')
    ax2.legend()
    
    # Add base labels on x-axis for zoomed plot (0-based indexing)
    ax2.set_xticks(zoom_positions)
    base_labels = []
    for pos in zoom_positions:
        if 0 <= pos < len(plasmid_seq):
            base = plasmid_seq[pos]  # 0-based indexing
        else:
            base = 'N'
        base_labels.append(f"{pos}\n{base}")
    ax2.set_xticklabels(base_labels, rotation=45, ha='right', fontsize=10)
    
    # Add coverage statistics for the zoomed region
    zoom_mean = np.mean(zoom_coverage)
    zoom_max = np.max(zoom_coverage)
    zoom_min = np.min(zoom_coverage)
    ax2.text(0.02, 0.98, f'Mean Coverage: {zoom_mean:.1f}\nMax Coverage: {zoom_max:.0f}\nMin Coverage: {zoom_min:.0f}', 
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the coverage plot
    coverage_plot_path = os.path.join(output_dir, "coverage_analysis.png")
    plt.savefig(coverage_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Coverage analysis plot saved: {coverage_plot_path}")
    
    # Save coverage data
    coverage_data_path = os.path.join(output_dir, "coverage_data.txt")
    with open(coverage_data_path, 'w') as fh:
        fh.write("Coverage Analysis\n")
        fh.write("=" * 50 + "\n")
        fh.write(f"Plasmid length: {plasmid_length} bp\n")
        fh.write(f"Total alignments: {len(summary_df)}\n")
        fh.write(f"Mean coverage: {mean_coverage:.2f}\n")
        fh.write(f"Max coverage: {max_coverage:.0f}\n")
        fh.write(f"Min coverage: {np.min(coverage):.0f}\n")
        fh.write(f"Standard deviation: {np.std(coverage):.2f}\n\n")
        
        fh.write(f"Sharpest Drop-off Analysis:\n")
        fh.write("-" * 30 + "\n")
        fh.write(f"Sharpest drop position: {zoom_center} (0-based)\n")
        fh.write(f"Drop magnitude: {drop_magnitude:.0f}\n")
        fh.write(f"Coverage at drop position: {coverage[zoom_center]:.0f}\n")
        fh.write(f"Coverage before drop: {coverage[max(0, zoom_center-1)]:.0f}\n")
        fh.write(f"Coverage after drop: {coverage[min(plasmid_length-1, zoom_center+1)]:.0f}\n")
        fh.write(f"Base at drop position: {plasmid_seq[zoom_center] if zoom_center < len(plasmid_seq) else 'N'}\n\n")
        
        fh.write(f"Zoomed Region (positions {zoom_start}-{zoom_end-1}, 0-based):\n")
        fh.write("-" * 40 + "\n")
        fh.write(f"Mean coverage: {zoom_mean:.2f}\n")
        fh.write(f"Max coverage: {zoom_max:.0f}\n")
        fh.write(f"Min coverage: {zoom_min:.0f}\n")
        fh.write(f"Standard deviation: {np.std(zoom_coverage):.2f}\n\n")
        
        fh.write("Position\tBase\tCoverage\n")
        fh.write("-" * 25 + "\n")
        for pos in zoom_positions:
            base = plasmid_seq[pos] if 0 <= pos < len(plasmid_seq) else 'N'
            fh.write(f"{pos}\t{base}\t{coverage[pos]:.0f}\n")
    
    logging.info(f"Coverage data saved: {coverage_data_path}")
    
    # Save coverage data as CSV for user replication
    import pandas as pd
    coverage_df = pd.DataFrame({
        'position': range(plasmid_length),
        'base': [plasmid_seq[i] if i < len(plasmid_seq) else 'N' for i in range(plasmid_length)],
        'coverage': coverage
    })
    
    coverage_csv_path = os.path.join(output_dir, "coverage_data.csv")
    coverage_df.to_csv(coverage_csv_path, index=False)
    logging.info(f"Coverage CSV saved: {coverage_csv_path}")
    
    # Save zoomed region data as CSV
    zoom_df = pd.DataFrame({
        'position': list(zoom_positions),
        'base': [plasmid_seq[pos] if 0 <= pos < len(plasmid_seq) else 'N' for pos in zoom_positions],
        'coverage': list(zoom_coverage)
    })
    
    zoom_csv_path = os.path.join(output_dir, "coverage_zoomed_data.csv")
    zoom_df.to_csv(zoom_csv_path, index=False)
    logging.info(f"Zoomed coverage CSV saved: {zoom_csv_path}")
    
    return coverage


def create_comprehensive_summary(summary_df, plasmid_length, output_dir):
    """Create comprehensive summary analysis across all samples with multiple testing correction."""
    
    # Overall position analysis
    start_counts = summary_df['start'].value_counts().sort_index()
    end_counts = summary_df['end'].value_counts().sort_index()
    
    total_alignments = len(summary_df)
    expected_frequency = total_alignments / plasmid_length
    
    # Statistical analysis with multiple testing correction
    start_significance = {}
    start_p_values = []
    start_positions = []
    
    for pos, count in start_counts.items():
        chi_square = (count - expected_frequency) ** 2 / expected_frequency
        p_value = 1 - scipy.stats.chi2.cdf(chi_square, 1)
        start_p_values.append(p_value)
        start_positions.append(pos)
        start_significance[pos] = {
            'count': count,
            'frequency': count / total_alignments,
            'expected': expected_frequency,
            'chi_square': chi_square,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
    
    # Apply multiple testing corrections for start positions
    if start_p_values:
        # Bonferroni correction
        bonferroni_threshold = 0.05 / len(start_p_values)
        
        for i, pos in enumerate(start_positions):
            # Bonferroni correction
            start_significance[pos]['bonferroni_p'] = min(start_p_values[i] * len(start_p_values), 1.0)
            start_significance[pos]['bonferroni_significant'] = start_p_values[i] < bonferroni_threshold
            
            # Benjamini-Hochberg correction (simplified)
            start_significance[pos]['bh_p'] = min(start_p_values[i] * len(start_p_values), 1.0)
            start_significance[pos]['bh_significant'] = start_p_values[i] < (0.05 * len(start_p_values) / (i + 1))
    
    end_significance = {}
    end_p_values = []
    end_positions = []
    
    for pos, count in end_counts.items():
        chi_square = (count - expected_frequency) ** 2 / expected_frequency
        p_value = 1 - scipy.stats.chi2.cdf(chi_square, 1)
        end_p_values.append(p_value)
        end_positions.append(pos)
        end_significance[pos] = {
            'count': count,
            'frequency': count / total_alignments,
            'expected': expected_frequency,
            'chi_square': chi_square,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
    
    # Apply multiple testing corrections for end positions
    if end_p_values:
        # Bonferroni correction
        bonferroni_threshold = 0.05 / len(end_p_values)
        
        for i, pos in enumerate(end_positions):
            # Bonferroni correction
            end_significance[pos]['bonferroni_p'] = min(end_p_values[i] * len(end_p_values), 1.0)
            end_significance[pos]['bonferroni_significant'] = end_p_values[i] < bonferroni_threshold
            
            # Benjamini-Hochberg correction (simplified)
            end_significance[pos]['bh_p'] = min(end_p_values[i] * len(end_p_values), 1.0)
            end_significance[pos]['bh_significant'] = end_p_values[i] < (0.05 * len(end_p_values) / (i + 1))
    
    # Create comprehensive summary plots with multiple testing correction
    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    
    # Plot 1: Overall start position distribution (Bonferroni)
    ax1 = axes[0, 0]
    positions = list(start_significance.keys())
    frequencies = [start_significance[pos]['frequency'] for pos in positions]
    bonferroni_sig = [start_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax1.bar(positions, frequencies, color=colors, alpha=0.7, edgecolor='black')
    ax1.axhline(y=1/plasmid_length, color='green', linestyle='--', alpha=0.7, label='Expected (uniform)')
    ax1.set_xlabel('Start Position (bp)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Overall Start Position Frequency\n(Red = Bonferroni Significant)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Overall end position distribution (Bonferroni)
    ax2 = axes[0, 1]
    positions = list(end_significance.keys())
    frequencies = [end_significance[pos]['frequency'] for pos in positions]
    bonferroni_sig = [end_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax2.bar(positions, frequencies, color=colors, alpha=0.7, edgecolor='black')
    ax2.axhline(y=1/plasmid_length, color='green', linestyle='--', alpha=0.7, label='Expected (uniform)')
    ax2.set_xlabel('End Position (bp)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Overall End Position Frequency\n(Red = Bonferroni Significant)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Start position frequency (BH)
    ax3 = axes[0, 2]
    positions = list(start_significance.keys())
    frequencies = [start_significance[pos]['frequency'] for pos in positions]
    bh_sig = [start_significance[pos]['bh_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bh_sig]
    ax3.bar(positions, frequencies, color=colors, alpha=0.7, edgecolor='black')
    ax3.axhline(y=1/plasmid_length, color='green', linestyle='--', alpha=0.7, label='Expected (uniform)')
    ax3.set_xlabel('Start Position (bp)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Overall Start Position Frequency\n(Red = BH Significant)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Start position Bonferroni p-values
    ax4 = axes[1, 0]
    positions = list(start_significance.keys())
    p_values = [start_significance[pos]['bonferroni_p'] for pos in positions]
    bonferroni_sig = [start_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax4.scatter(positions, p_values, color=colors, alpha=0.7, s=20)
    ax4.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    ax4.set_xlabel('Start Position (bp)')
    ax4.set_ylabel('Bonferroni P-value')
    ax4.set_title('Start Position Bonferroni P-values')
    ax4.set_yscale('log')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: End position Bonferroni p-values
    ax5 = axes[1, 1]
    positions = list(end_significance.keys())
    p_values = [end_significance[pos]['bonferroni_p'] for pos in positions]
    bonferroni_sig = [end_significance[pos]['bonferroni_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bonferroni_sig]
    ax5.scatter(positions, p_values, color=colors, alpha=0.7, s=20)
    ax5.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    ax5.set_xlabel('End Position (bp)')
    ax5.set_ylabel('Bonferroni P-value')
    ax5.set_title('End Position Bonferroni P-values')
    ax5.set_yscale('log')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Start position BH p-values
    ax6 = axes[1, 2]
    positions = list(start_significance.keys())
    p_values = [start_significance[pos]['bh_p'] for pos in positions]
    bh_sig = [start_significance[pos]['bh_significant'] for pos in positions]
    
    colors = ['red' if sig else 'blue' for sig in bh_sig]
    ax6.scatter(positions, p_values, color=colors, alpha=0.7, s=20)
    ax6.axhline(y=0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    ax6.set_xlabel('Start Position (bp)')
    ax6.set_ylabel('BH P-value')
    ax6.set_title('Start Position BH P-values')
    ax6.set_yscale('log')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax4.set_yscale('log')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Alignment length distribution
    ax5 = axes[1, 1]
    ax5.hist(summary_df['length'], bins=50, alpha=0.7, edgecolor='black', color='green')
    ax5.set_xlabel('Alignment Length (bp)')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Overall Alignment Length Distribution')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Mapping quality distribution
    ax6 = axes[1, 2]
    ax6.hist(summary_df['mapping_quality'], bins=30, alpha=0.7, edgecolor='black', color='red')
    ax6.set_xlabel('Mapping Quality')
    ax6.set_ylabel('Frequency')
    ax6.set_title('Overall Mapping Quality Distribution')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    summary_plot_path = os.path.join(output_dir, "comprehensive_summary_plot.png")
    plt.savefig(summary_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Comprehensive summary plot saved: {summary_plot_path}")
    
    # Save comprehensive analysis with multiple testing corrections
    summary_analysis_path = os.path.join(output_dir, "comprehensive_position_analysis.txt")
    with open(summary_analysis_path, 'w') as fh:
        fh.write("COMPREHENSIVE POSITION ANALYSIS\n")
        fh.write("=" * 60 + "\n\n")
        
        # Multiple testing correction information
        fh.write("MULTIPLE TESTING CORRECTION\n")
        fh.write("-" * 35 + "\n")
        fh.write(f"Number of start positions tested: {len(start_significance)}\n")
        fh.write(f"Number of end positions tested: {len(end_significance)}\n")
        fh.write(f"Bonferroni threshold (α=0.05): {0.05/len(start_significance):.6f}\n")
        fh.write(f"Benjamini-Hochberg FDR control: α=0.05\n\n")
        
        fh.write("OVERALL STATISTICS\n")
        fh.write("-" * 20 + "\n")
        fh.write(f"Total alignments across all samples: {total_alignments}\n")
        fh.write(f"Plasmid length: {plasmid_length} bp\n")
        fh.write(f"Expected frequency (uniform): {1/plasmid_length:.6f}\n")
        fh.write(f"Average alignment length: {summary_df['length'].mean():.2f} bp\n")
        fh.write(f"Average mapping quality: {summary_df['mapping_quality'].mean():.2f}\n\n")
        
        # Most significant start positions
        sorted_starts = sorted(start_significance.items(), key=lambda x: x[1]['count'], reverse=True)
        fh.write("MOST COMMON START POSITIONS (All Samples)\n")
        fh.write("-" * 45 + "\n")
        fh.write("Position\tCount\tFrequency\tRaw P-value\tBonferroni P\tBH P\tBonferroni Sig\tBH Sig\n")
        fh.write("-" * 100 + "\n")
        
        for pos, stats in sorted_starts[:30]:  # Top 30
            fh.write(f"{pos}\t{stats['count']}\t{stats['frequency']:.6f}\t{stats['p_value']:.6f}\t{stats['bonferroni_p']:.6f}\t{stats['bh_p']:.6f}\t{stats['bonferroni_significant']}\t{stats['bh_significant']}\n")
        
        fh.write("\n\nMOST COMMON END POSITIONS (All Samples)\n")
        fh.write("-" * 43 + "\n")
        sorted_ends = sorted(end_significance.items(), key=lambda x: x[1]['count'], reverse=True)
        fh.write("Position\tCount\tFrequency\tRaw P-value\tBonferroni P\tBH P\tBonferroni Sig\tBH Sig\n")
        fh.write("-" * 100 + "\n")
        
        for pos, stats in sorted_ends[:30]:  # Top 30
            fh.write(f"{pos}\t{stats['count']}\t{stats['frequency']:.6f}\t{stats['p_value']:.6f}\t{stats['bonferroni_p']:.6f}\t{stats['bh_p']:.6f}\t{stats['bonferroni_significant']}\t{stats['bh_significant']}\n")
        
        # Statistical summary with multiple testing corrections
        bonferroni_sig_starts = sum(1 for stats in start_significance.values() if stats['bonferroni_significant'])
        bonferroni_sig_ends = sum(1 for stats in end_significance.values() if stats['bonferroni_significant'])
        bh_sig_starts = sum(1 for stats in start_significance.values() if stats['bh_significant'])
        bh_sig_ends = sum(1 for stats in end_significance.values() if stats['bh_significant'])
        
        fh.write(f"\n\nSTATISTICAL SUMMARY (Multiple Testing Corrected)\n")
        fh.write("-" * 45 + "\n")
        fh.write(f"Bonferroni significant start positions: {bonferroni_sig_starts}\n")
        fh.write(f"Bonferroni significant end positions: {bonferroni_sig_ends}\n")
        fh.write(f"BH significant start positions: {bh_sig_starts}\n")
        fh.write(f"BH significant end positions: {bh_sig_ends}\n")
        fh.write(f"Total unique start positions: {len(start_significance)}\n")
        fh.write(f"Total unique end positions: {len(end_significance)}\n")
    
    logging.info(f"Comprehensive analysis saved: {summary_analysis_path}")
    
    # Create coverage analysis (merged paired-end data)
    logging.info("Creating coverage analysis with merged paired-end data...")
    create_coverage_analysis(summary_df, plasmid_length, output_dir)
    
    logging.info("Plasmid mapping analysis completed successfully!")

if __name__ == "__main__":
    main() 