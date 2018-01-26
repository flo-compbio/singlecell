#!/usr/bin/env python3

"""Script to map inDrop reads using STAR.

"""

import sys
import argparse

from genometools import misc

from .. import mapping

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = 'Map inDrop reads using STAR.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Input parameters')

    g.add_argument(
        '-r', '--read-file', type=str, required=True,
        help='FASTQ file (uncompressed) containing the reads.')

    g.add_argument(
        '-i', '--index-dir', type=str, required=True,
        help='Directory containing the STAR index.')

    g = parser.add_argument_group('Output parameters')        

    g.add_argument(
        '-os', '--output-script-file', type=str, required=True,
        help='Output file containing the bash script used for running STAR.')

    g.add_argument(
        '-ol', '--output-log-file', type=str, required=True,
        help='Text file containing the output messages produced by STAR.')

    g.add_argument(
        '-od', '--output-mapping-dir', type=str, required=True,
        help='STAR output directory.')

    g = parser.add_argument_group('Other parameters')        

    g.add_argument(
        '-d', '--use-docker', action='store_true',
        help='Whether to use a dockerized version of STAR.'
    )

    return parser


def main(args=None):
    """Entry point for script."""
    
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    read_file = args.read_file
    index_dir = args.index_dir
    output_script_file = args.output_script_file
    output_log_file = args.output_log_file
    output_mapping_dir = args.output_mapping_dir
    use_docker = args.use_docker

    misc.make_sure_dir_exists(output_mapping_dir)

    mapping.map_with_star(read_file, index_dir,
                          output_script_file, output_log_file,
                          output_mapping_dir,
                          compressed=False, use_docker=use_docker)


    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
