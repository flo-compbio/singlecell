#!/usr/bin/env python3

"""Create a configuration file template for use with the inDrop pipeline."""

import sys
import argparse

from genometools import misc

from .. import config

_LOGGER = misc.get_logger()


def get_argument_parser():

    desc = ('Create a standard configuration file (in YAML) format for use '
            'with the inDrop pipeline.')

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Output parameters')

    g.add_argument(
        '-o', '--output-file', type=str, required=True,
        help='The output file.')    

    return parser


def main(args=None):
    """Entry point for script."""
    
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    output_file = args.output_file

    config_str = config.get_config_template()

    with open(output_file, 'w') as ofh:
        ofh.write(config_str)

    _LOGGER.info('Written configuration file template to "%s".', output_file)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
