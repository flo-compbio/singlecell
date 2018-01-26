#!/usr/bin/env python3

"""Script for counting transcriptomic reads for each inDrop barcode."""

import argparse

from genometools import misc

from .. import barcodes

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = 'Count transcriptomic reads for each inDrop barcode.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-a', '--alignment-file', type=str, required=True,
        help='.bam file containing the mapped reads.')

    g.add_argument(
        '-b1', '--barcode1-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-b2', '--barcode2-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-g', '--gene-file', type=str, required=True,
        help='.tsv file containing list of protein-coding genes.')

    g.add_argument(
        '-n', '--genome-annotation-file', type=str, required=True,
        help='.gtf file containing genome annotations.')

    g.add_argument(
        '-o', '--output-file', type=str, required=True,
        help='Output file.')

    g = parser.add_argument_group('Counting parameters')

    g.add_argument(
        '-m', '--max-genes', type=int, required=False, default=0,
        help=('Maxmimum number of genes to process. If 0, '
              'process all genes. [0]')
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
    gene_file = args.gene_file
    genome_annotation_file = args.genome_annotation_file
    output_file = args.output_file
    max_genes = args.max_genes

    if max_genes == 0:
        max_genes = None

    barcodes.count_transcriptomic_reads(
        alignment_file, barcode1_file, barcode2_file,
        gene_file, genome_annotation_file,
        output_file, max_genes=max_genes)

