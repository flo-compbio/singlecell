#!/usr/bin/env python3

"""Script for counting mapped reads for each inDrop barcode."""

import sys
import time
import argparse

import pandas as pd
import numpy as np
import pysam

from genometools import misc

from .. import barcodes

_LOGGER = misc.get_logger()


def get_argument_parser():
    desc = 'Count mapped reads for all inDrop barcodes.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-a', '--alignment-file', type=str, required=True,
        help='.bam file containing the aligned reads.')

    g.add_argument(
        '-b1', '--barcode1-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-b2', '--barcode2-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-o', '--output-file', type=str, required=True,
        help='Output file.')

    g = parser.add_argument_group('Counting parameters')

    g.add_argument(
        '-m', '--max-reads', type=int, required=False, default=0,
        help='Maxmimum number of reads to process. If 0, '
             'process all reads. [0]'
    )

    return parser


def main(args=None):
    """Entry point for script."""
    
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    alignment_file = args.alignment_file
    barcode1_file = args.barcode1_file
    barcode2_file = args.barcode2_file
    output_file = args.output_file

    max_reads = args.max_reads
    if max_reads == 0:
        max_reads = None

    barcodes.count_mapped_reads(
        alignment_file, barcode1_file, barcode2_file,
        output_file,
        max_reads=max_reads)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)