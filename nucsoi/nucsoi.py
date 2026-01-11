#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
from tqdm import tqdm

__version__ = "1.0.0"

def setup_directories(output_dir):
    """Create necessary directories for the NuCSOI pipeline."""
    dirs = {
        'inputs_raw': os.path.join(output_dir, 'inputs', 'raw_fastqgzs'),
        'inputs_qc': os.path.join(output_dir, 'inputs', 'qc_reads'),
        'inputs_plasmid': os.path.join(output_dir, 'inputs', 'plasmid'),
        'results': os.path.join(output_dir, 'results'),
        'scripts': os.path.join(output_dir, 'scripts')
    }
    
    print("Creating directory structure...")
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
    
    return dirs

def copy_scripts(scripts_dir):
    """Copy the necessary scripts to the output directory."""
    # Get the package directory
    package_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Copy required scripts
    scripts_to_copy = [
        'scripts/map_to_plasmid.py',
        'scripts/fastp.py'
    ]
    
    print("Copying analysis scripts...")
    for script in tqdm(scripts_to_copy, desc="Copying scripts"):
        src = os.path.join(package_dir, script)
        dst = os.path.join(scripts_dir, os.path.basename(script))
        if os.path.exists(src):
            shutil.copy2(src, dst)
        else:
            print(f"Warning: {src} not found")

def create_configs_yaml(output_dir, quality_cutoff):
    """Create configs.yaml with the specified quality cutoff."""
    config_content = f"""# NuCSOI Configuration File
# ========================

# Input/Output Directories
input_dir_raw: inputs/raw_fastqgzs
output_dir_qc_reads: inputs/qc_reads
results_base: results

# Quality Control Parameters
quality_threshold: {quality_cutoff}

# Plasmid Mapping Parameters
# (These are automatically detected from inputs/plasmid/)
"""
    
    config_path = os.path.join(output_dir, 'configs.yaml')
    with open(config_path, 'w') as f:
        f.write(config_content)
    
    print(f"Created configs.yaml with quality threshold: {quality_cutoff}")

def create_makefile(output_dir):
    """Create Makefile for the NuCSOI pipeline."""
    makefile_content = """# NuCSOI Pipeline Makefile
# ======================

.PHONY: all clean qc_reads map_plasmid

all: qc_reads map_plasmid

qc_reads:
	python scripts/fastp.py -c configs.yaml

map_plasmid:
	python scripts/map_to_plasmid.py -c configs.yaml

clean:
	rm -rf results/*
	rm -rf inputs/qc_reads/*
	rm -rf inputs/paired/*
	find . -name "*.log" -delete
	find . -name "*.tmp" -delete
"""
    
    makefile_path = os.path.join(output_dir, 'Makefile')
    with open(makefile_path, 'w') as f:
        f.write(makefile_content)
    
    print("Created Makefile")

def copy_fastq_files(fastq_files, raw_dir):
    """Copy FASTQ files to the raw directory."""
    print("Copying FASTQ files...")
    for fastq_file in tqdm(fastq_files, desc="Copying FASTQ files"):
        if not os.path.exists(fastq_file):
            print(f"Error: FASTQ file not found: {fastq_file}")
            return False
        
        dst = os.path.join(raw_dir, os.path.basename(fastq_file))
        shutil.copy2(fastq_file, dst)
    
    return True

def copy_plasmid_reference(plasmid_ref, plasmid_dir):
    """Copy plasmid reference to the plasmid directory."""
    if not os.path.exists(plasmid_ref):
        print(f"Error: Plasmid reference not found: {plasmid_ref}")
        return False
    
    print("Copying plasmid reference...")
    dst = os.path.join(plasmid_dir, os.path.basename(plasmid_ref))
    shutil.copy2(plasmid_ref, dst)
    print(f"Copied {plasmid_ref} to {dst}")
    
    return True

def run_pipeline(output_dir):
    """Run the complete NuCSOI pipeline."""
    print("\n" + "="*60)
    print("Running NuCSOI Pipeline")
    print("="*60)
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        # Run the complete pipeline with real-time output
        print("Executing pipeline: make all")
        print("Progress:")
        
        # Run make with real-time output
        process = subprocess.Popen(
            ['make', 'all'],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # Read output in real-time
        for line in iter(process.stdout.readline, ''):
            if line.strip():
                # Filter and format the output
                if 'INFO:' in line or 'Processing' in line or 'Running' in line:
                    # Extract the important part of the log message
                    if 'INFO:' in line:
                        # Extract the message after INFO:
                        message = line.split('INFO:')[-1].strip()
                        if message:
                            print(f"  {message}")
                    else:
                        print(f"  {line.strip()}")
        
        process.stdout.close()
        return_code = process.wait()
        
        if return_code == 0:
            print("\nNuCSOI pipeline completed successfully!")
            print(f"\nResults are available in: {output_dir}")
            print("\nKey output files:")
            print("- results/plasmid_mapping_*/coverage_analysis.png")
            print("- results/plasmid_mapping_*/coverage_data.txt")
            print("- results/plasmid_mapping_*/comprehensive_summary_plot.png")
            print("- results/plasmid_mapping_*/comprehensive_position_analysis.txt")
        else:
            print("\nNuCSOI pipeline failed!")
            return False
            
    except Exception as e:
        print(f"Error running pipeline: {e}")
        return False
    finally:
        os.chdir(original_dir)
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description="NuCSOI: Nuclease Cleavage Site and Overhang Identification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  nucsoi -f sample_R1.fastq.gz sample_R2.fastq.gz -p plasmid.fasta -o results/
  nucsoi -f sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz -p plasmid.fasta -o results/ -q 30
  nucsoi -f sample_R1.fastq.gz sample_R2.fastq.gz -p plasmid.fasta -o results/ --run-pipeline

For more information, visit: https://github.com/Matt115A/NuCSOI
        """
    )
    
    parser.add_argument(
        '-f', '--fastq-files',
        nargs='+',
        required=True,
        help='Paired FASTQ files (R1 and R2 files)'
    )
    
    parser.add_argument(
        '-p', '--plasmid',
        required=True,
        help='Plasmid reference FASTA file'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '-q', '--quality-cutoff',
        type=int,
        default=30,
        help='Quality cutoff for fastp (default: 30)'
    )
    
    parser.add_argument(
        '--run-pipeline',
        action='store_true',
        help='Automatically run the pipeline after setup'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'NuCSOI {__version__}'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if len(args.fastq_files) % 2 != 0:
        print("Error: Number of FASTQ files must be even (paired R1/R2 files)")
        sys.exit(1)
    
    # Create output directory
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    print("NuCSOI Pipeline Setup")
    print("===================")
    print(f"Output directory: {output_dir}")
    print(f"Quality cutoff: {args.quality_cutoff}")
    print(f"FASTQ files: {len(args.fastq_files)} files")
    print(f"Plasmid reference: {args.plasmid}")
    print()
    
    # Setup directories
    dirs = setup_directories(output_dir)
    
    # Copy scripts
    copy_scripts(dirs['scripts'])
    
    # Create configuration files
    create_configs_yaml(output_dir, args.quality_cutoff)
    create_makefile(output_dir)
    
    # Copy input files
    if not copy_fastq_files(args.fastq_files, dirs['inputs_raw']):
        sys.exit(1)
    
    if not copy_plasmid_reference(args.plasmid, dirs['inputs_plasmid']):
        sys.exit(1)
    
    print(f"\nSetup complete! Files copied to: {output_dir}")
    
    if args.run_pipeline:
        if not run_pipeline(output_dir):
            sys.exit(1)
    else:
        print(f"\nTo run the pipeline, navigate to {output_dir} and run:")
        print("  make all")
        print("\nOr run with --run-pipeline flag to execute automatically.")

if __name__ == "__main__":
    main() 