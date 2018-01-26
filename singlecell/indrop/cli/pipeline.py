#!/usr/bin/env python3

"""inDrop pipeline script.

"""

import sys
import argparse

from genometools import misc

from .. import pipeline

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = 'Process inDrop reads.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Parameters')

    g.add_argument(
        '-c', '--config-file', type=str, required=True,
        help='Pipeline configuration file (in YAML format).')

    return parser


def main(args=None):
    """Entry point for script."""
    
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    config_file = args.config_file

    pipeline.run_pipeline(config_file)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
