# NuCSI Pip Installation Guide

## Overview

This guide explains how to make NuCSI installable via pip and what files are needed.

## Required Files for Pip Installation

### 1. Package Structure
```
NUclease-Cleavage-Site-Id/
├── nucsi/                          # Main package directory
│   ├── __init__.py                 # Package initialization
│   ├── nucsi.py                    # Main script
│   └── scripts/                    # Analysis scripts
│       ├── fastp.py
│       └── map_to_plasmid.py
├── setup.py                        # Traditional setup script
├── pyproject.toml                  # Modern Python packaging
├── MANIFEST.in                     # Package manifest
├── LICENSE                         # MIT license
├── README.md                       # Documentation
└── requirements.txt                # Dependencies
```

### 2. Key Files Created/Modified

#### `setup.py`
- Traditional Python packaging setup
- Defines package metadata, dependencies, and entry points
- Handles package discovery and data files

#### `pyproject.toml`
- Modern Python packaging configuration
- Defines build system requirements
- Specifies project metadata and dependencies
- Handles entry points and package data

#### `MANIFEST.in`
- Specifies which files to include in the package
- Includes README, LICENSE, and scripts
- Excludes cache files

#### `nucsi/__init__.py`
- Package initialization file
- Defines version and imports
- Exports main function

#### `LICENSE`
- MIT license for the project
- Required for PyPI distribution

## Installation Methods

### 1. Local Development Installation
```bash
# Clone the repository
git clone https://github.com/Matt115A/NuCSI.git
cd NuCSI

# Install in development mode
pip install -e .
```

### 2. Local Package Installation
```bash
# Build the package
python setup.py sdist bdist_wheel

# Install from wheel
pip install dist/nucsi-1.0.0-py3-none-any.whl
```

### 3. PyPI Installation (Future)
```bash
# Once published to PyPI
pip install nucsi
```

## Usage After Installation

### Command Line Interface
```bash
# Basic usage
nucsi -f sample_R1.fastq.gz sample_R2.fastq.gz -p plasmid.fasta -o results/

# With quality cutoff
nucsi -f sample_R1.fastq.gz sample_R2.fastq.gz -p plasmid.fasta -o results/ -q 30

# Automatic pipeline execution
nucsi -f sample_R1.fastq.gz sample_R2.fastq.gz -p plasmid.fasta -o results/ --run-pipeline
```

### Python Import
```python
import nucsi
from nucsi import main

# Access version
print(nucsi.__version__)
```

## Publishing to PyPI

### 1. Prerequisites
- PyPI account (https://pypi.org)
- Test PyPI account (https://test.pypi.org)
- `twine` package installed

### 2. Build and Upload
```bash
# Build the package
python setup.py sdist bdist_wheel

# Upload to Test PyPI first
twine upload --repository testpypi dist/*

# Test installation from Test PyPI
pip install --index-url https://test.pypi.org/simple/ nucsi

# Upload to PyPI
twine upload dist/*
```

### 3. Version Management
- Update version in `nucsi/__init__.py`
- Tag releases in Git
- Use semantic versioning (e.g., 1.0.0, 1.0.1, 1.1.0)

## Dependencies

### Required Python Packages
- pandas>=1.3.0
- numpy>=1.21.0
- scipy>=1.7.0
- pyyaml>=5.4.0
- tqdm>=4.62.0
- biopython>=1.79.0
- pysam>=0.19.0
- matplotlib>=3.4.0
- seaborn>=0.11.0
- logomaker>=0.8.0

### External Tools (Not Included)
- fastp (quality control)
- bwa (read mapping)
- samtools (BAM processing)

## Troubleshooting

### Common Issues

1. **Package not found**: Ensure all files are in the correct locations
2. **Import errors**: Check `__init__.py` and package structure
3. **Scripts not found**: Verify `MANIFEST.in` includes scripts
4. **Dependencies missing**: Check `requirements.txt` and `setup.py`

### Testing Installation
```bash
# Test package import
python -c "import nucsi; print('Success')"

# Test command line tool
nucsi --help

# Test version
nucsi --version
```

## Next Steps

1. **Test thoroughly** with different Python versions
2. **Add CI/CD** for automated testing and deployment
3. **Create documentation** for PyPI
4. **Set up automated releases** with GitHub Actions
5. **Monitor package usage** and issues

## Notes

- The package uses both `setup.py` and `pyproject.toml` for maximum compatibility
- Scripts are included as package data for easy access
- Entry points create the `nucsi` command-line tool
- All dependencies are specified in both `requirements.txt` and `pyproject.toml` 