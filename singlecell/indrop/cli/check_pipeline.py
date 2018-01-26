#!/usr/bin/env python3

"""inDrop pipeline checking script.

"""

import os
import sys
import argparse
import subprocess

import yaml

from genometools import misc

from ... import util
from .. import config

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


def check_if_process_ok(cmd):
    subproc = subprocess.Popen(cmd, shell=True)
    subproc.communicate()
    return subproc.returncode == 0


def file_exists(path):
    exists = True
    try:
        with open(path) as fh:
            pass
    except FileNotFoundError:
        exists = False
    return exists


def file_is_readable(path):
    is_readable = True
    try:
        with open(path) as fh:
            pass
    except IOError:
        is_readable = False
    return is_readable


def main(args=None):
    """Entry point for script."""
    
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    config_file = args.config_file

    if not os.path.isfile(config_file):
        _LOGGER.error('Config file not found! Specify a config file and then '
                      're-run this script.')
        return 1

    valid = True

    conf, errors = config.read_config(config_file)
    if errors:
        _LOGGER.info('There were errors loading the configuration file.')
        valid = False
    else:
        _LOGGER.info('Successfully loaded the configuration file.')

    input_ = conf['input']
    output = conf['output']
    params = conf['parameters']
    pipeline = conf['pipeline']
    output_dir = output['output_dir']

    #barcode1_file = resource_filename(
    #    'singlecell', 'data/indrop/gel_barcode1_list.txt')
    #barcode2_file = resource_filename(
    #    'singlecell', 'data/indrop/gel_barcode2_list.txt-')

    # check if output directory is valid
    if not os.path.isdir(output_dir):
        parent_dir = os.path.dirname(os.path.abspath(output_dir))
        if not os.path.isdir(parent_dir):
            _LOGGER.error('[output/output_dir]: '
                          'The output directory "%s" does not exist, and its '
                          'parent directory does not exist either. You need to'
                          ' make sure that the parent directory of the '
                          'specified output director exists!', output_dir)
            valid = False
        else:
            _LOGGER.info('The output directory does not exist, but its parent '
                         'directory does, so the pipeline will automatically '
                         'create the output directory when it runs.')
    else:
        if not util.is_empty_dir(output_dir):
            if not output['allow_nonempty_output_dir']:
                _LOGGER.error('[output/output_dir], '
                              '[output/allow_nonempty_output_dir]: '
                              'The output directory is not empty, but '
                              '"allow_nonempty_output_dir: yes" has not been '
                              'specified. You need '
                              'to either make sure the output directory is '
                              'empty, or specify '
                              '"allow_nonempty_output_dir: yes".')
                valid = False
            else:
                _LOGGER.info('The output directory is not empty, but '
                             '"allow_nonempty_output_dir: yes" has been '
                             'specified in the configuration file, so that\'s '
                             'OK.')
        else:
            _LOGGER.info('The output directory already exists, but it\'s '
                         'empty, so that works.')

    # check if STAR and samtools are executable, or if docker is present,
    # if applicable
    if params['use_docker']:
        # check if docker is working
        if not check_if_process_ok('docker run hello-world'):
            _LOGGER.error('[parameters/use_docker]: '
                          'You have configured the pipeline to run using '
                          'docker, but there was a problem with running '
                          '"docker run hello-world". Make sure docker is set '
                          'up correctly!')
            valid = False
        else:
            _LOGGER.info('docker seems to be installed.')
    else:
        # check if STAR and samtools are present
        if not check_if_process_ok('STAR --version'):
            _LOGGER.error('There was a problem with running "STAR --version". '
                          'Make sure the STAR executable is in your PATH. '
                          'Alternatively, you might want to consider running '
                          'the pipeline using docker (parameters/use_docker).')
            valid = False
        else:
            _LOGGER.info('STAR seems to be installed.')
        if not check_if_process_ok('samtools --version'):
            _LOGGER.error('There was a problem with running '
                          '"samtools --version". '
                          'Make sure the samtools executable is in your PATH.'
                          'Alternatively, you might want to consider running '
                          'the pipeline using docker (parameters/use_docker).')
            valid = False
        else:
            _LOGGER.info('samtools seems to be installed.')

    def check_file(path, option, name):
        valid_file = True
        if not file_exists(path):
            _LOGGER.error('[%s]: '
                          'The %s file "%s" does not exist. Make sure you '
                          'specify the correct file path in the configuration '
                          'file.', option, name, path)
            valid_file = False

        elif not file_is_readable(path):
            _LOGGER.error('The %s file "%s" exists, but is not readable. Make '
                          'sure you have permission to read the file')
            valid_file = False
        else:
            _LOGGER.info('The %s file is readable!', name)
        return valid_file

    if not check_file(
            input_['barcode_read_file'],
            'input/barcode_read_file',
            'barcode read (FASTQ.gz)'):
        valid = False

    if not check_file(
            input_['mrna_read_file'],
            'input/mrna_read_file',
            'mRNA read (FASTQ.gz)'):
        valid = False

    if not check_file(
            input_['genome_file'],
            'input/genome_file',
            'genome sequence (.fa.gz)'):
        valid = False

    if not check_file(
            input_['genome_annotation_file'],
            'input/genome_annotation_file',
            'genome annotation (.gtf.gz)'):
        valid = False

    if not os.path.isdir(input_['star_index_dir']):
        _LOGGER.error('[input/star_index_dir]: '
                      'The STAR index directory does not exist! Make sure you '
                      'specify the correct directory path in the configuration'
                      ' file.')
        valid = False
    elif not os.access(input_['star_index_dir'], os.R_OK | os.X_OK):
        _LOGGER.error('[input/star_index_dir]: '
                      'The STAR index directory exists, but is not readable! '
                      'Make sure you have permission to read the directory '
                      'contents.')
        valid = False

    _LOGGER.info('')
    if not valid:
        _LOGGER.warning('There have been problems identified with the setup / '
                        'configuration of your inDrop pipeline. Fix these '
                        'problems before running the pipeline, in order to '
                        'avoid a pipeline crash.')
    else:
        _LOGGER.info('No problems found!')

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
