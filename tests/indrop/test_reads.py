"""Test functions for the `reads` module."""

# Author: Florian Wagner <florian.wagner@nyu.edu>

import os
import logging

from pkg_resources import resource_filename
import pytest

from genometools import misc

from singlecell import indrop

_LOGGER = logging.getLogger(__name__)


@pytest.mark.online
def test_process_reads(my_barcode_read_file, my_mrna_read_file,
                       my_output_pypath):
    """Tests the read processing function."""
    
    barcode1_file = resource_filename(
        'singlecell', 'data/indrop/gel_barcode1_list.txt')
    barcode2_file = resource_filename(
        'singlecell', 'data/indrop/gel_barcode2_list.txt')

    output_read_file = str(my_output_pypath.join('processed_reads.fastq'))
    output_count_file = str(my_output_pypath.join('barcode_counts.tsv'))
    output_log_file = str(my_output_pypath.join('read_processing_log.txt'))
    
    indrop.reads.process_reads(
        my_barcode_read_file, my_mrna_read_file,
        barcode1_file, barcode2_file,
        output_read_file, output_count_file, output_log_file)

    md5sum = misc.get_file_md5sum(output_read_file)
    assert md5sum == '4a21825d92ed776f126dd89843ba16d6'

    md5sum = misc.get_file_md5sum(output_count_file)
    assert md5sum == 'b2776081a49442bb01e62e1470574f18'
