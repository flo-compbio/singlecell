#!/usr/bin/env python3

"""Generate a STAR index for use with the inDrop pipeline.

"""

import sys
import os
import argparse
import logging
import time

from genometools import misc
from genometools import ensembl

from .. import mapping

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = 'Generate a STAR index for use with the inDrop pipeline.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Input parameters')

    g.add_argument(
        '-g', '--decompressed-genome-file', type=str, required=True,
        help='.fasta file containing the genome sequence.')

    g.add_argument(
        '-n', '--decompressed-genome-annotation-file', type=str, required=True,
        help='.gtf file containing the genome annotations.')

    g = parser.add_argument_group('Output parameters')

    g.add_argument(
        '-od', '--output-dir', type=str, required=True,
        help='Output directory for the index.')    

    g.add_argument(
        '-os', '--output-script-file', type=str, required=True,
        help=('Output file containing the command executed to generate the '
              'index.'))

    g.add_argument(
        '-ol', '--output-log-file', type=str, required=True,
        help='Output file containing the index.')    

    g = parser.add_argument_group('Other parameters')

    g.add_argument(
        '-t', '--num-threads', type=int, required=False, default=16,
        help='Number of threads to use. [16]')

    g.add_argument(
        '-d', '--use-docker', action='store_true',
        help='Whether to use a dockerized version of STAR.'
    )

    return parser


def main(args=None):
    """Entry point for script."""
    
    t0 = time.time()

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    decompressed_genome_file = args.decompressed_genome_file
    decompressed_genome_annotation_file = \
            args.decompressed_genome_annotation_file

    output_dir = args.output_dir
    output_script_file = args.output_script_file
    output_log_file = args.output_log_file

    num_threads = args.num_threads
    use_docker = args.use_docker

    misc.make_sure_dir_exists(output_dir)
    
    mapping.generate_star_index(
        decompressed_genome_file, decompressed_genome_annotation_file,
        output_dir, output_script_file, output_log_file,
        num_threads=num_threads, use_docker=use_docker)

    #_LOGGER.removeHandler(file_handler)
    t1 = time.time()
    t = t1 - t0
    _LOGGER.info('Index generation finished in %.1f s (%.1f min)!', t, t/60)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
