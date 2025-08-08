# NuCSI: Nuclease Cleavage Site Identification Tool

A command-line interface for identifying nuclease cleavage sites from paired-end sequencing data with single base-pair resolution.

## Overview

NuCSI is a bioinformatics pipeline that processes paired-end sequencing data to identify nuclease cleavage sites on plasmid references. The tool performs quality control, plasmid mapping, and comprehensive coverage analysis to provide statistical insights into nuclease activity with single base-pair resolution.

The pipeline is specifically designed for analyzing nuclease cleavage experiments where:
- DNA is cleaved by a nuclease at specific sites
- The cleaved DNA is sequenced using paired-end sequencing
- Coverage analysis reveals cleavage sites as positions with reduced read coverage
- Statistical analysis identifies significant cleavage patterns

## Features

- **Quality Control**: Fastp-based quality filtering with configurable thresholds
- **Plasmid Mapping**: BWA-based mapping to circular plasmid references
- **Coverage Analysis**: Automatic detection of coverage drop-offs and zoomed analysis
- **Statistical Analysis**: Multiple hypothesis testing correction (Bonferroni and Benjamini-Hochberg)
- **Position Analysis**: Single base-pair resolution analysis of mapping start/end positions
- **Data Reproducibility**: CSV exports for all coverage and position data
- **Visualization**: Comprehensive plots and coverage analysis
- **Circular Plasmid Support**: Handles circular plasmid references with proper wrapping

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
- Should contain sequencing data from nuclease-treated samples

### Plasmid Reference
- Must be in FASTA format
- Should be the circular plasmid reference sequence
- Will be automatically normalized to uppercase
- Should be the same plasmid used in the nuclease experiment

## Output Structure

```
output_dir/
├── inputs/
│   ├── raw_fastqgzs/          # Copied input FASTQ files
│   ├── qc_reads/              # Quality-controlled reads
│   └── plasmid/               # Plasmid reference files
├── results/
│   └── plasmid_mapping_*/     # Results for each plasmid
│       ├── coverage_analysis.png
│       ├── coverage_data.csv          # Full coverage data (reproducible)
│       ├── coverage_zoomed_data.csv   # Zoomed region data (reproducible)
│       ├── coverage_data.txt
│       ├── comprehensive_summary_plot.png
│       └── comprehensive_position_analysis.txt
├── scripts/                   # Analysis scripts
├── configs.yaml              # Configuration file
└── Makefile                  # Pipeline makefile
```

## Key Output Files

### Coverage Analysis
- `coverage_analysis.png`: Coverage plots for entire plasmid and zoomed regions around cleavage sites
- `coverage_data.csv`: Full coverage data for all positions (reproducible)
- `coverage_zoomed_data.csv`: Zoomed region data around sharpest drop-off (reproducible)
- `coverage_data.txt`: Detailed coverage statistics and base information

### Position Analysis
- `comprehensive_summary_plot.png`: Statistical analysis plots with multiple testing correction
- `comprehensive_position_analysis.txt`: Detailed position statistics with Bonferroni and Benjamini-Hochberg corrections

### Data Reproducibility
All coverage and position data are exported as CSV files, allowing users to:
- Recreate plots using their preferred visualization tools
- Perform additional statistical analyses
- Integrate data into other pipelines
- Share reproducible datasets

## Pipeline Steps

1. **Quality Control**: Fastp processing with specified quality cutoff (default: Q30)
2. **Plasmid Mapping**: Map reads to plasmid reference using BWA with circular plasmid support
3. **Coverage Analysis**: Calculate and visualize coverage with automatic drop-off detection
4. **Position Analysis**: Statistical analysis of mapping positions with multiple testing correction
5. **Cleavage Site Detection**: Identify positions with lowest coverage as potential cleavage sites

## Biological Interpretation

### Coverage Patterns
- **High coverage regions**: Intact DNA with many reads mapping
- **Low coverage regions**: Potential cleavage sites with reduced read mapping
- **Sharp coverage drops**: Indicate nuclease cleavage activity

### Cleavage Site Identification
- **Primary method**: Position with lowest coverage
- **Secondary method**: Position with sharpest coverage drop
- **Statistical significance**: Multiple testing correction applied
- **Base context**: Sequence context provided for each cleavage site

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
- `fastp`: Quality control and adapter trimming
- `bwa`: Read mapping to reference sequences
- `samtools`: BAM file processing and indexing

### Python Packages
- `pandas`: Data manipulation and analysis
- `numpy`: Numerical computing
- `scipy`: Statistical analysis
- `matplotlib`: Plotting and visualization
- `seaborn`: Statistical data visualization
- `biopython`: Bioinformatics tools
- `pysam`: SAM/BAM file processing
- `tqdm`: Progress bars
- `pyyaml`: YAML configuration file processing

## Troubleshooting

### Common Issues

1. **FASTQ files not found**: Ensure files exist and paths are correct
2. **Plasmid reference not found**: Check file path and format
3. **Pipeline fails**: Check dependencies are installed (fastp, bwa, samtools)
4. **Memory issues**: Reduce quality cutoff or use smaller datasets
5. **No alignments found**: Check plasmid reference matches your sequencing data

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

