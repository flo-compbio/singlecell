#!/usr/bin/env python3

"""Script to process raw inDrop reads and count barcodes.

"""

import sys
import argparse

from genometools import misc

from .. import reads

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = 'Process inDrop reads.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-rb', '--barcode-read-file', type=str, required=True,
        help='.fastq.gz file containing the reads with barcode sequences.')

    g.add_argument(
        '-rm', '--mrna-read-file', type=str, required=True,
        help='.fastq.gz file containing the reads with mRNA sequences.')

    g.add_argument(
        '-b1', '--barcode1-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-b2', '--barcode2-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-or', '--output-read-file', type=str, required=True,
        help='Output read (FASTQ) file.')

    g.add_argument(
        '-oc', '--output-count-file', type=str, required=True,
        help='Output count file.')

    g.add_argument(
        '-ol', '--output-log-file', type=str, required=True,
        help='Output log file.')

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
    
    barcode_read_file = args.barcode_read_file
    mrna_read_file = args.mrna_read_file
    barcode1_file = args.barcode1_file
    barcode2_file = args.barcode2_file
    output_read_file = args.output_read_file
    output_count_file = args.output_count_file
    output_log_file = args.output_log_file

    max_reads = args.max_reads
    if max_reads == 0:
        max_reads = None

    reads.process_reads(
        barcode_read_file, mrna_read_file,
        barcode1_file, barcode2_file,
        output_read_file, output_count_file, output_log_file,
        max_reads=max_reads)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
