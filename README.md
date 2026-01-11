# NuCSOI

Nuclease Cleavage Site and Overhang Identification

## Description

NuCSOI processes paired-end sequencing data to identify nuclease cleavage sites on plasmid references. The pipeline performs quality control, mapping, coverage analysis, and statistical analysis to identify cleavage sites with single base-pair resolution.

## Installation

```bash
pip install nucsoi
```

External dependencies (install separately):
- fastp (quality control)
- bwa (read mapping)
- samtools (BAM processing)

## Inputs

- **FASTQ files**: Paired-end sequencing files (R1 and R2). Must be an even number of files.
- **Plasmid reference**: FASTA file containing the circular plasmid sequence.

## Usage

```bash
nucsoi -f R1.fastq.gz R2.fastq.gz -p plasmid.fasta -o output_dir/
```

Use `nucsoi --help` for all available options.

## Pipeline Stages

1. **Quality Control**: Filters reads using fastp with configurable quality threshold (default: Q30).
2. **Plasmid Mapping**: Maps quality-filtered reads to the plasmid reference using BWA. Handles circular references.
3. **Coverage Analysis**: Calculates coverage at each position. Identifies regions with coverage drop-offs.
4. **Position Analysis**: Statistical analysis of mapping positions. Applies multiple testing correction (Bonferroni and Benjamini-Hochberg).

## Outputs

Results are written to the specified output directory:

```
output_dir/
├── inputs/
│   ├── raw_fastqgzs/          # Input FASTQ files
│   ├── qc_reads/              # Quality-controlled reads
│   └── plasmid/               # Plasmid reference
├── results/
│   └── plasmid_mapping_*/
│       ├── coverage_analysis.png
│       ├── coverage_data.csv
│       ├── coverage_zoomed_data.csv
│       ├── coverage_data.txt
│       ├── comprehensive_summary_plot.png
│       └── comprehensive_position_analysis.txt
├── scripts/                   # Analysis scripts
├── configs.yaml              # Configuration file
└── Makefile                  # Pipeline makefile
```

- `coverage_analysis.png`: Coverage plots for entire plasmid and zoomed regions
- `coverage_data.csv`: Coverage data for all positions
- `coverage_zoomed_data.csv`: Coverage data for zoomed regions
- `coverage_data.txt`: Coverage statistics
- `comprehensive_summary_plot.png`: Statistical analysis plots
- `comprehensive_position_analysis.txt`: Position statistics with multiple testing corrections

## Options

- `-f, --fastq-files`: Paired FASTQ files (required)
- `-p, --plasmid`: Plasmid reference FASTA file (required)
- `-o, --output-dir`: Output directory (required)
- `-q, --quality-cutoff`: Quality cutoff for fastp (default: 30)
- `--run-pipeline`: Automatically run pipeline after setup
- `--version`: Show version number
- `-h, --help`: Show help message

## License

MIT License
