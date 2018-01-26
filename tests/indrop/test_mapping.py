"""Test functions for the `indrop.reads` module."""

import logging

# import pytest

from singlecell.indrop import mapping

_LOGGER = logging.getLogger(__name__)


def test_generate_STAR_index(my_dummy_dir, my_dummy_file):
    """Tests the STAR index generation function.
    
    We're not actually running STAR, because it requires so much memory. We're
    only calling the index generation function and tell it to skip running
    STAR."""
    
    #output_dir = str(tmpdir_factory.mktemp('singlecell_mapping_output',
    #                                       numbered=False))
    
    # make sure the index_dir is an existing directory
    output_dir = str(my_dummy_dir)

    decompressed_genome_file = my_dummy_file
    decompressed_genome_annotation_file = my_dummy_file

    output_script_file = str(my_dummy_dir.join('map_with_star.sh'))
    output_log_file = str(my_dummy_dir.join('mapping_log.txt'))

    mapping.generate_star_index(
        decompressed_genome_file, decompressed_genome_annotation_file,
        output_dir, output_script_file, output_log_file,
        num_threads=1, use_docker=False, simulate=True)


def test_mapping_with_STAR(my_dummy_dir, my_dummy_file):
    """Tests the STAR mapping function.
    
    We're not actually running STAR, because it requires so much memory. We're
    only calling the mapping function and tell it to skip running STAR."""
    
    # make sure the index_dir is an existing directory
    index_dir = str(my_dummy_dir)
    read_file = my_dummy_file

    output_script_file = str(my_dummy_dir.join('map_with_star.sh'))
    output_log_file = str(my_dummy_dir.join('mapping_log.txt'))
    output_mapping_dir = str(my_dummy_dir.join('mapping'))

    mapping.map_with_star(
        read_file, index_dir,
        output_script_file, output_log_file, output_mapping_dir,
        num_threads=1, compressed=False, use_docker=True,
        simulate=True)
