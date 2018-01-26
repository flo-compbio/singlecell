"""Test functions for the `barcodes` module."""

import logging

from pkg_resources import resource_filename
import pytest

from genometools import misc

from singlecell.indrop import barcodes

_LOGGER = logging.getLogger(__name__)

@pytest.mark.online
def test_count_mapped(my_alignment_file, my_output_pypath):
    """Tests the read processing function."""

    barcode1_file = resource_filename(
        'singlecell', 'data/indrop/gel_barcode1_list.txt')
    barcode2_file = resource_filename(
        'singlecell', 'data/indrop/gel_barcode2_list.txt')

    output_file = str(my_output_pypath.join('barcode_counts_mapped.tsv'))
    barcodes.count_mapped_reads(
        my_alignment_file, barcode1_file, barcode2_file, output_file)
    
    md5sum = misc.get_file_md5sum(output_file)
    assert md5sum == '9b137a6ee925f02273cefb08bd9e1366'
