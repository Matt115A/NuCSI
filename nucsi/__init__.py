#!/usr/bin/env python3
"""
NuCSI: Nuclease Cleavage Site Identification Tool

A command-line interface for identifying nuclease cleavage sites from paired-end sequencing data.
"""

__version__ = "1.0.0"
__author__ = "Matthew Penner"
__email__ = "matthew.penner@crick.ac.uk"

from .nucsi import main

__all__ = ["main"] 