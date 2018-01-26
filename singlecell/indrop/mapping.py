"""Functions for mapping inDrop reads."""

import os
import subprocess

from jinja2 import Environment, PackageLoader, select_autoescape

from genometools import misc

from .. import util

_LOGGER = misc.get_logger()

_TEMPLATE_ENV = Environment(
    loader=PackageLoader('singlecell',
                         os.path.join('data', 'templates')),
    autoescape=select_autoescape(['html', 'xml'])
)


def generate_star_index(
        decompressed_genome_file, decompressed_genome_annotation_file,
        output_dir, output_script_file, output_log_file,
        num_threads=1, use_docker=False, simulate=False):
    """Generates a STAR index.
    
    TODO: docstring"""

    assert os.path.isfile(decompressed_genome_file)
    assert os.path.isfile(decompressed_genome_annotation_file)
    
    misc.make_sure_dir_exists(output_dir)
    template = _TEMPLATE_ENV.get_template(
        os.path.join('star', 'build_index.sh'))

    script = template.render(
        decompressed_genome_file=decompressed_genome_file,
        decompressed_genome_annoation_file=decompressed_genome_annotation_file,
        output_dir=output_dir,
        num_threads=num_threads,
        use_docker=use_docker)
    
    with open(output_script_file, 'w') as ofh:
        ofh.write(script)
    util.make_file_executable(output_script_file)

    output_script_file = os.path.abspath(output_script_file)

    _LOGGER.info('Starting to build a STAR index...')

    if not simulate:
        subproc = subprocess.Popen('%s 2>&1 > %s'
                                % (output_script_file, output_log_file),
                                shell=True)
        subproc.communicate()
        _LOGGER.info('STAR index build process exited with return code %d',
                    subproc.returncode)
    
    _LOGGER.info('Done building the STAR index!')   


def map_with_star(read_file, index_dir,
                  output_script_file, output_log_file,
                  output_mapping_dir,
                  num_threads=1, compressed=True, use_docker=False,
                  simulate=False,
                  **kwargs):
    """Maps the mRNA sequence-containing read to the genome using STAR.
    
    TODO: docstring"""
    
    assert os.path.isfile(read_file)
    assert os.path.isdir(index_dir)

    misc.make_sure_dir_exists(output_mapping_dir)
    
    output_mapping_prefix = output_mapping_dir.rstrip(os.sep) + os.sep

    template = _TEMPLATE_ENV.get_template(
        os.path.join('star', 'map_single-end.sh'))

    script = template.render(
        read_file=read_file,
        index_dir=index_dir,
        output_prefix=output_mapping_prefix,
        compressed=compressed,
        out_sam_type='BAM SortedByCoordinate',
        out_mult_nmax=1,
        num_threads=num_threads,
        use_docker=use_docker,
        **kwargs
    )
    
    with open(output_script_file, 'w') as ofh:
        ofh.write(script)
        
    util.make_file_executable(output_script_file)

    output_script_file = os.path.abspath(output_script_file)

    _LOGGER.info('Starting to map with STAR...')
    if not simulate:
        subproc = subprocess.Popen('%s 2>&1 > %s'
                                   % (output_script_file, output_log_file),
                                   shell=True)
        subproc.communicate()
        _LOGGER.info('Mapping process exited with return code %d',
                     subproc.returncode)

    _LOGGER.info('Done with mapping!')
