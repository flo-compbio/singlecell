#cython: language_level=3, boundscheck=False

"""Functions for working with and analyzing inDrop barcodes."""

import logging
import sys
import time

import pandas as pd
import numpy as np
import pysam

_LOGGER = logging.getLogger(__name__)

from .. import util


def count_mapped_reads(alignment_file, barcode1_file, barcode2_file,
                       output_file,
                       max_reads=None):
    """Count the number of mapped reads for each inDrop barcode.
    
    TODO: docstring"""

    # matrix for storing numbers for all barcode combinations
    cdef unsigned long int[:, ::1] counts = np.zeros((384, 384), dtype=np.uint64)

    # first barcodes
    barcodes1 = pd.read_csv(
        barcode1_file, sep='\t', header=None, squeeze=True) \
            .apply(util.get_reverse_complement).values
    assert barcodes1.ndim == 1
    
    # second barcodes
    barcodes2 = pd.read_csv(
        barcode2_file, sep='\t', header=None, squeeze=True) \
            .apply(util.get_reverse_complement).values
    assert barcodes2.ndim == 1
    
    cdef int bc1, bc2
    cdef str query_name
    cdef unsigned long int c_max_reads, i
    
    if max_reads is None:
        c_max_reads = 0
    else:
        c_max_reads = max_reads
    
    t0 = time.time()
    with pysam.AlignmentFile(alignment_file, 'rb') as bam:
        # assert bam.check_index()
        i = 0
        for read in bam.fetch():
            if c_max_reads != 0 and i >= c_max_reads:
                break
            i += 1
            query_name = read.query_name
            bc1 = int(query_name[:3])
            bc2 = int(query_name[4:7])
            counts[bc1, bc2] += 1

            #if (i+1) % 1000000 == 0:
            #    print(i+1, end=' ')
            #    sys.stdout.flush()

    t1 = time.time()
    print()
    t = t1-t0
    _LOGGER.info('Processed %d alignments in %.1f s.', i, t)
    _LOGGER.info("(That's %.1f s per million reads.)", 1e6*(t/i))

    bc1_index = pd.Index(barcodes1)
    bc2_index = pd.Index(barcodes2)
    df = pd.DataFrame(np.uint64(counts), index=bc1_index, columns=bc2_index)
    df.to_csv(output_file, sep='\t')
