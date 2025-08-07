# NuCSI: Nuclease Cleavage Site Identification Tool

A command-line interface for identifying nuclease cleavage sites from paired-end sequencing data.

## Overview

NuCSI is a bioinformatics pipeline that processes paired-end sequencing data to identify nuclease cleavage sites on plasmid references. The tool performs quality control, sequence scanning, plasmid mapping, and coverage analysis to provide comprehensive insights into nuclease activity.

## Features

- **Quality Control**: Fastp-based quality filtering with configurable thresholds
- **Sequence Scanning**: Automatic detection of sequences between specified motifs
- **Plasmid Mapping**: BWA-based mapping to circular plasmid references
- **Coverage Analysis**: Automatic detection of coverage drop-offs and zoomed analysis
- **Statistical Analysis**: Multiple hypothesis testing correction and position analysis
- **Visualization**: Comprehensive plots and coverage analysis

## Quick Start

### Installation

1. Clone the repository:
```bash
git clone https://github.com/Matt115A/NuCSI.git
cd NuCSI
```

2. Install dependencies:
```bash
# Create conda environment
conda env create -f environment.yml
conda activate nucsi

# Or install manually
pip install -r requirements.txt
conda install -c bioconda fastp bwa samtools
```

### Basic Usage

```bash
# Single sample analysis
nucsi.py -f sample_R1.fastq.gz sample_R2.fastq.gz -p plasmid.fasta -o results/

# Multiple samples with automatic execution
nucsi.py -f sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz \
         -p plasmid.fasta -o results/ -q 30 --run-pipeline
```

## Command Line Options

| Option | Description | Required | Default |
|--------|-------------|----------|---------|
| `-f, --fastq-files` | Paired FASTQ files (R1 and R2 files) | Yes | - |
| `-p, --plasmid` | Plasmid reference FASTA file | Yes | - |
| `-o, --output-dir` | Output directory for results | Yes | - |
| `-q, --quality-cutoff` | Quality cutoff for fastp | No | 30 |
| `--run-pipeline` | Automatically run the pipeline after setup | No | False |
| `--version` | Show program's version number | No | - |
| `-h, --help` | Show help message | No | - |

## Input Requirements

### FASTQ Files
- Must be paired-end sequencing data (R1 and R2 files)
- Files should be in `.fastq.gz` format
- Number of files must be even (paired R1/R2 files)

### Plasmid Reference
- Must be in FASTA format
- Should be the circular plasmid reference sequence
- Will be automatically normalized to uppercase

## Output Structure

```
output_dir/
├── inputs/
│   ├── raw_fastqgzs/          # Copied input FASTQ files
│   ├── qc_reads/              # Quality-controlled reads
│   └── plasmid/               # Plasmid reference files
├── results/
│   ├── plasmid_mapping_*/     # Results for each plasmid
│   │   ├── coverage_analysis.png
│   │   ├── coverage_data.txt
│   │   ├── comprehensive_summary_plot.png
│   │   └── comprehensive_position_analysis.txt
│   └── consensus/             # Sequence consensus results
├── scripts/                   # Analysis scripts
├── configs.yaml              # Configuration file
└── Makefile                  # Pipeline makefile
```

## Key Output Files

### Coverage Analysis
- `coverage_analysis.png`: Coverage plots for entire plasmid and zoomed regions
- `coverage_data.txt`: Detailed coverage statistics and base information

### Position Analysis
- `comprehensive_summary_plot.png`: Statistical analysis plots
- `comprehensive_position_analysis.txt`: Detailed position statistics

### Sequence Analysis
- `consensus/`: Multiple sequence alignment and logo plots

## Pipeline Steps

1. **Quality Control**: Fastp processing with specified quality cutoff
2. **Sequence Scanning**: Identify sequences between specified motifs
3. **Plasmid Mapping**: Map reads to plasmid reference
4. **Coverage Analysis**: Calculate and visualize coverage
5. **Position Analysis**: Statistical analysis of mapping positions

## Examples

### Single Sample Analysis
```bash
nucsi.py -f sample_R1.fastq.gz sample_R2.fastq.gz \
         -p plasmid.fasta \
         -o results/ \
         -q 30 \
         --run-pipeline
```

### Multiple Sample Analysis
```bash
nucsi.py -f sample1_R1.fastq.gz sample1_R2.fastq.gz \
         sample2_R1.fastq.gz sample2_R2.fastq.gz \
         -p plasmid.fasta \
         -o results/ \
         -q 25 \
         --run-pipeline
```

### Setup Only (Manual Execution)
```bash
nucsi.py -f sample_R1.fastq.gz sample_R2.fastq.gz \
         -p plasmid.fasta \
         -o results/ \
         -q 30

cd results/
make all
```

## Dependencies

### Required Tools
- `fastp`: Quality control
- `bwa`: Read mapping
- `samtools`: BAM file processing
- `mafft`: Multiple sequence alignment

### Python Packages
- `pandas`, `numpy`, `scipy`
- `matplotlib`, `seaborn`
- `biopython`, `pysam`
- `tqdm`, `pyyaml`

## Troubleshooting

### Common Issues

1. **FASTQ files not found**: Ensure files exist and paths are correct
2. **Plasmid reference not found**: Check file path and format
3. **Pipeline fails**: Check dependencies are installed (fastp, bwa, samtools)
4. **Memory issues**: Reduce quality cutoff or use smaller datasets

### Getting Help

1. Check the troubleshooting section
2. Review the GitHub issues
3. Create a new issue with detailed information

## Citation

If you use NuCSI in your research, please cite:

```
NuCSI: Nuclease Cleavage Site Identification Tool
https://github.com/Matt115A/NuCSI
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

