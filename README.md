# NuCSI (Nuclease Cleavage Site Identification)

A comprehensive bioinformatics pipeline for identifying and analyzing nuclease cleavage sites from high-throughput sequencing data.

## Overview

NuCSI is designed to process Illumina sequencing data and identify nuclease cleavage sites through multiple analysis approaches:

1. **Quality Control**: Process raw FASTQ files using fastp
2. **Sequence Scanning**: Identify sequences between specified motifs
3. **Plasmid Mapping**: Map reads to reference plasmids and analyze alignment positions

## Features

- **Quality Control**: Automated FASTQ processing with configurable quality thresholds
- **Motif-Based Sequence Extraction**: Extract sequences between upstream and downstream motifs
- **Plasmid Mapping**: Map reads to circular plasmid references with statistical analysis
- **Position Analysis**: Single-base resolution analysis of alignment start/end positions
- **Multiple Testing Correction**: Bonferroni and Benjamini-Hochberg corrections for statistical significance
- **Visualization**: Comprehensive plotting including zoom plots for top hit regions

## Installation

### Prerequisites

- Python 3.10+
- Conda package manager
- Unix-like operating system (Linux/macOS)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/Matt115A/NuCSI.git
   cd NuCSI
   ```

2. Run the setup script:
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```

3. Activate the conda environment:
   ```bash
   conda activate nucsi
   ```

## Usage

### Quick Start

Run the complete pipeline:
```bash
make all
```

This will execute:
1. Quality control of raw FASTQ files
2. Sequence scanning between motifs
3. Plasmid mapping and position analysis

### Individual Steps

#### 1. Quality Control
```bash
make qc_reads
# or
python scripts/fastp.py -c configs.yaml
```

Processes raw FASTQ files in `inputs/raw_fastqgzs/` and outputs quality-controlled reads to `inputs/qc_reads/`.

#### 2. Sequence Scanning
```bash
make scan_sequences
# or
python scripts/scan_sequences.py -c configs.yaml
```

Scans quality-controlled reads for sequences between specified upstream and downstream motifs. Results are saved in `results/consensus/`.

#### 3. Plasmid Mapping
```bash
make map_plasmid
# or
python scripts/map_to_plasmid.py -c configs.yaml
```

Maps reads to reference plasmids and performs detailed position analysis. Results are saved in `results/plasmid_mapping/`.

### Configuration

Edit `configs.yaml` to customize the analysis:

```yaml
# Input/Output Directories
input_dir_raw: inputs/raw_fastqgzs
output_dir_qc_reads: inputs/qc_reads
results_base: results

# Quality Control Parameters
quality_threshold: 30

# Sequence Scanning Parameters
upstream: ccagcttcaaaa
downstream: gcgctctgaagtt
```

## Input Requirements

### Directory Structure
```
NuCSI/
├── inputs/
│   ├── raw_fastqgzs/          # Raw FASTQ files
│   ├── qc_reads/              # Quality-controlled reads (auto-generated)
│   └── plasmid/               # Reference plasmid FASTA files
├── results/                   # Analysis results (auto-generated)
└── scripts/                   # Analysis scripts
```

### File Naming Conventions

- **Raw FASTQ files**: `{sample}_R1_001.fastq.gz` and `{sample}_R2_001.fastq.gz`
- **Plasmid references**: Any `.fasta` file in `inputs/plasmid/`

## Output Files

### Quality Control
- Quality-controlled FASTQ files in `inputs/qc_reads/`
- Fastp reports in `results/logs/`

### Sequence Scanning
- Consensus sequences and alignments in `results/consensus/`
- Logo plots showing sequence conservation

### Plasmid Mapping
- BAM files and alignment statistics in `results/plasmid_mapping/`
- Position analysis plots and tables
- Zoom plots for top hit regions (20 bases on either side)
- Comprehensive summary analysis

## Key Outputs

### Position Analysis
- **Top Start Position**: Position 255 (366 alignments, 27.6%)
- **Top End Position**: Position 251 (396 alignments, 29.9%)
- **Statistical Significance**: Bonferroni and Benjamini-Hochberg corrected p-values
- **Zoom Plots**: Single-base resolution visualization of top hit regions

### Sequence Analysis
- Extracted sequences between specified motifs
- Multiple sequence alignments
- Consensus logos and conservation plots

## Statistical Analysis

The pipeline includes comprehensive statistical analysis:

- **Chi-square tests** against uniform distribution
- **Multiple testing correction** (Bonferroni and Benjamini-Hochberg)
- **Single-base resolution** position analysis
- **Significance thresholds** with FDR control

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Ensure conda environment is activated
2. **File permissions**: Make sure setup.sh is executable
3. **Memory issues**: Large files may require increased memory allocation

### Logs

All analysis steps generate detailed logs in `results/logs/` with timestamps.

## Citation

If you use NuCSI in your research, please cite:

```
NuCSI: Nuclease Cleavage Site Identification Pipeline
[Your citation information here]
```

## License

[Add your license information here]

## Contributing

[Add contribution guidelines here]

