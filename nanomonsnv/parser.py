#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .detect import detect_main
from .add_control import add_control_main
from .filter import filter_main

def create_parser():
    prog = "nanomonsnv"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
    
    def _create_detect_parser(subparsers):
        
        detect_parser = subparsers.add_parser("detect", help = "Detecting somatic mutation candidates using Fisher exact tests")
        detect_parser.add_argument("tumor_bam", type = str, help = "Path to tumor bam file")
        detect_parser.add_argument("control_bam", type = str, help = "Path to matched control file")
        detect_parser.add_argument("output_file", type = str, help = "Path to the output file")
        detect_parser.add_argument("reference", type = str, help = "Path to the reference genome")
        return detect_parser

    def _create_add_control_parser(subparsers):

        add_control_parser = subparsers.add_parser("add_control", help = "Adding non-matched control infomation to the variant candidate file")
        add_control_parser.add_argument("variant_file", type = str, help = "Path to variant candidate file")
        add_control_parser.add_argument("output_file", type = str, help = "Path to the output file")
        add_control_parser.add_argument("reference", type = str, help = "Path to the reference genome")
        add_control_parser.add_argument("control_bams", type = str, nargs='+', help = "Path to non-matched control files")
        return add_control_parser
        
    def _create_filt_parser(subparsers):
    
        filt_parser = subparsers.add_parser("filter", help = "Validation using Illumina short read tumor and control sequencing data")
        filt_parser.add_argument("variant_file", type = str, help = "Path to variant candidate file")
        filt_parser.add_argument("tumor_bam", type = str, help = "Path to Illumina tumor bam file")
        filt_parser.add_argument("control_bam", type = str, help = "Path to Illumina matched control file")
        filt_parser.add_argument("output_file", type = str, help = "Path to the output file")
        filt_parser.add_argument("reference", type = str, help = "Path to the reference genome")
        return filt_parser
        
    detect_parser = _create_detect_parser(subparsers)
    detect_parser.set_defaults(func = detect_main)
    add_control_parser = _create_add_control_parser(subparsers)
    add_control_parser.set_defaults(func = add_control_main)
    filt_parser = _create_filt_parser(subparsers)
    filt_parser.set_defaults(func = filter_main)
    return parser
