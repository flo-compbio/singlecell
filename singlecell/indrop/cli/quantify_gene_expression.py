#!/usr/bin/env python3

"""Gene expression quantification script for inDrop data."""

import sys
import argparse

from genometools import misc

from .. import expression

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = 'Quantify gene expression based on aligned inDrop reads.'

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
        '-c', '--barcode-count-file', type=str, required=True,
        help='.tsv file containing the read counts for all barcodes.')

    g.add_argument(
        '-g', '--gene-file', type=str, required=True,
        help='.tsv file containing list of protein-coding genes.')

    g.add_argument(
        '-n', '--genome-annotation-file', type=str, required=True,
        help='.gtf file containing genome annotations.')

    g.add_argument(
        '-or', '--raw-output-file', type=str, required=True,
        help='Output file containing the raw (unfiltered) gene read counts.')

    g.add_argument(
        '-of', '--filtered-output-file', type=str, required=True,
        help='Output file containing the UMI-filtered gene read counts.')

    g = parser.add_argument_group('Other arguments')

    g.add_argument(
        '-p', '--cell-prefix', type=str, default='',
        help='Prefix for cell names. (default: no prefix)')
    
    g = parser.add_argument_group('Counting parameters')

    g.add_argument(
        '--num-cells', type=int, required=True,
        help=('Number of cells (barcodes) included in the analysis.')
    )
    
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
    barcode_count_file = args.barcode_count_file
    gene_file = args.gene_file
    genome_annotation_file = args.genome_annotation_file
    raw_output_file = args.raw_output_file
    filtered_output_file = args.filtered_output_file
    num_cells = args.num_cells
    max_genes = args.max_genes
    cell_prefix = args.cell_prefix

    if max_genes == 0:
        max_genes = None

    expression.quantify_gene_expression(
        alignment_file, barcode_count_file,
        gene_file, genome_annotation_file,
        raw_output_file, filtered_output_file,
        num_cells,
        cell_prefix=cell_prefix, max_genes=max_genes)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
