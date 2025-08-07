# NuCSI: Nuclease Cleavage Site Identification Tool

A command-line interface for identifying nuclease cleavage sites from paired-end sequencing data.

## Overview

NuCSI is a bioinformatics pipeline that processes paired-end sequencing data to identify nuclease cleavage sites on plasmid references. The tool performs quality control, plasmid mapping, and comprehensive coverage analysis to provide statistical insights into nuclease activity with single base-pair resolution.

## Features

- **Quality Control**: Fastp-based quality filtering with configurable thresholds
- **Plasmid Mapping**: BWA-based mapping to circular plasmid references
- **Coverage Analysis**: Automatic detection of coverage drop-offs and zoomed analysis
- **Statistical Analysis**: Multiple hypothesis testing correction (Bonferroni and Benjamini-Hochberg)
- **Position Analysis**: Single base-pair resolution analysis of mapping start/end positions
- **Data Reproducibility**: CSV exports for all coverage and position data
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
- `coverage_analysis.png`: Coverage plots for entire plasmid and zoomed regions
- `coverage_data.csv`: Full coverage data for all positions (reproducible)
- `coverage_zoomed_data.csv`: Zoomed region data around sharpest drop-off (reproducible)
- `coverage_data.txt`: Detailed coverage statistics and base information

### Position Analysis
- `comprehensive_summary_plot.png`: Statistical analysis plots
- `comprehensive_position_analysis.txt`: Detailed position statistics with multiple testing correction

### Data Reproducibility
All coverage and position data are exported as CSV files, allowing users to:
- Recreate plots using their preferred visualization tools
- Perform additional statistical analyses
- Integrate data into other pipelines
- Share reproducible datasets

## Pipeline Steps

1. **Quality Control**: Fastp processing with specified quality cutoff
2. **Plasmid Mapping**: Map reads to plasmid reference using BWA
3. **Coverage Analysis**: Calculate and visualize coverage with automatic drop-off detection
4. **Position Analysis**: Statistical analysis of mapping positions with multiple testing correction

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

## Accessibility and Distribution

### Current Availability
- **GitHub Repository**: Available at https://github.com/Matt115A/NuCSI
- **Installation**: Manual installation via git clone and conda/pip
- **Documentation**: Comprehensive README and inline code documentation

### Enhanced Accessibility Options

#### 1. **Conda Package Distribution**
```bash
# Future: Install via conda
conda install -c bioconda nucsi
```

#### 2. **PyPI Package Distribution**
```bash
# Future: Install via pip
pip install nucsi
```

#### 3. **Docker Container**
```dockerfile
# Future: Dockerfile for containerized execution
FROM continuumio/miniconda3
RUN conda install -c bioconda fastp bwa samtools
RUN pip install pandas numpy scipy matplotlib seaborn biopython pysam tqdm pyyaml
COPY nucsi.py /usr/local/bin/
RUN chmod +x /usr/local/bin/nucsi.py
ENTRYPOINT ["nucsi.py"]
```

#### 4. **Galaxy Tool Integration**
- Integration with Galaxy workflow system
- Web-based interface for non-command-line users
- Workflow sharing and reproducibility

#### 5. **Web Application**
- Flask/FastAPI-based web interface
- File upload and processing
- Real-time progress tracking
- Result visualization

## Dependencies

### Required Tools
- `fastp`: Quality control
- `bwa`: Read mapping
- `samtools`: BAM file processing

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

