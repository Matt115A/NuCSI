#!/usr/bin/env python3
"""
Setup script for NuCSOI: Nuclease Cleavage Site and Overhang Identification
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        requirements = []
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                requirements.append(line)
        return requirements

# Get version from nucsoi package
def get_version():
    with open("nucsoi/__init__.py", "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip('"').strip("'")
    return "1.0.0"

setup(
    name="nucsoi",
    version=get_version(),
    author="Matthew Penner",
    author_email="matthew.penner@crick.ac.uk",
    description="Nuclease Cleavage Site and Overhang Identification",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/Matt115A/NuCSOI",
    project_urls={
        "Bug Reports": "https://github.com/Matt115A/NuCSOI/issues",
        "Source": "https://github.com/Matt115A/NuCSOI",
        "Documentation": "https://github.com/Matt115A/NuCSOI#readme",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10,<3.11",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "nucsoi": ["scripts/*.py"],
    },
    install_requires=read_requirements(),
    entry_points={
        "console_scripts": [
            "nucsoi=nucsoi.nucsoi:main",
        ],
    },
    keywords="bioinformatics, nuclease, cleavage, sequencing, coverage, analysis",
    zip_safe=False,
) 