#cython: language_level=3, boundscheck=False

"""Functions for working with inDrop reads."""

#from libc.stdio  cimport FILE, fopen, fread, fclose, fgets, sscanf
#from libc.string cimport strlen, strcmp
from libc.stdio  cimport FILE, fopen, fclose, fgets, fputs
from libc.string cimport strlen
from libc.stdlib cimport malloc

#mport 
cimport numpy as np
np.import_array()

import logging
import time
import io
import tempfile
import subprocess
import os

import pandas as pd
import numpy as np

from singlecell import util
from singlecell.indrop import barcodes

_LOGGER = logging.getLogger(__name__)


def create_gzip_pipe(path):

    tfh = tempfile.NamedTemporaryFile(mode='wb', prefix='indrop_',
                                      delete=False)
    temp = tfh.name
    _LOGGER.debug('Temp file: %s', temp)
    tfh.close()
    os.remove(temp) # delete the temp file
    os.mkfifo(temp)
    subproc = subprocess.Popen('gunzip -c "%s" > "%s"' %(path, temp),
                               shell=True)
    return temp



def process_reads(barcode_read_file, mrna_read_file,
                  barcode1_file, barcode2_file,
                  output_read_file, output_count_file,
                  max_reads=None):
    """Process inDrop sequencing reads.

    The inDrop method relies on paired-end sequencing. One mate contains
    information about the cell identity (defined by a specific pair of
    barcodes), as well as the transcript identity (defined by a specific UMI
    sequence). The other mate contains the 3'-end sequence of an mRNA molecule
    (just before the poly-A tail).

    This function tries to identify the barcode pair of each read, and counts
    the number of reads with each barcode combination in the process. For reads
    where the barcode can be successfully identified, it adds the barcode
    information, as well as the UMI index, to the read name of the mRNA mate.
    In this way, the other mate is then no longer needed in downstream steps of
    the pipeline.

    TODO: docstring
    """

    w1 = 'GAGTGATTGCTTGTGACGCCTT'
    
    w1_mapping = {w1}
    w1_mapping |= set(util.get_edit_sequences(w1, 1))
    w1_mapping |= set(util.get_edit_sequences(w1, 2))

    bc1_exact_mapping = barcodes.get_exact_mapping(barcode1_file)
    bc1_mismatch_mapping = barcodes.get_mismatch_mapping(barcode1_file)
    
    bc2_exact_mapping = barcodes.get_exact_mapping(barcode2_file)
    bc2_mismatch_mapping = barcodes.get_mismatch_mapping(barcode2_file)   
    
    # first barcodes
    barcodes1 = pd.read_csv(
        barcode1_file, sep='\t', header=None, squeeze=True) \
            .apply(util.get_reverse_complement).values
    assert barcodes1.ndim == 1
    cdef unsigned long long int[::1] bc1_counts = np.zeros(barcodes1.size, dtype=np.uint64)
    
    # second barcodes
    barcodes2 = pd.read_csv(
        barcode2_file, sep='\t', header=None, squeeze=True) \
            .apply(util.get_reverse_complement).values
    assert barcodes2.ndim == 1
    cdef unsigned long long int[::1] bc2_counts = np.zeros(barcodes2.size, dtype=np.uint64)
    
    cdef unsigned long long int[:, ::1] bc_counts = np.zeros((barcodes1.size, barcodes2.size), dtype=np.uint64)
        
    pipe1 = create_gzip_pipe(barcode_read_file)
    pipe2 = create_gzip_pipe(mrna_read_file)
    _LOGGER.info('Named pipe 1: %s', pipe1)
    _LOGGER.info('Named pipe 2: %s', pipe2)
    
    cdef char* ret
    cdef int buf_size = 1000
    # cdef size_t nl = <size_t>newline_chars
    cdef char* buf = <char*>malloc(buf_size)
    cdef char* seq = <char*>malloc(buf_size)
    cdef char* qual = <char*>malloc(buf_size)
    # cdef char* l1 = <char*>malloc(buf_size)
    cdef char* l2 = <char*>malloc(buf_size)
    cdef char* l3 = <char*>malloc(buf_size)
    cdef char* l4 = <char*>malloc(buf_size)
    cdef FILE* file1
    cdef FILE* file2
    cdef int bc1, bc2, bc1_len, r1_len, good_qual
    cdef int bc1_found, bc2_found, bc1_start, w1_len, bc2_start
    cdef int p, umi_start
    cdef unsigned long long int passed_w1, passed_bc, passed_qual
    cdef unsigned long long int c, c_max_reads
    
    if max_reads is None:
        c_max_reads = 0
    else:
        c_max_reads = max_reads
    
    pipe1_bytes = pipe1.encode('ascii')
    file1 = fopen(pipe1_bytes, 'r')
    
    pipe2_bytes = pipe2.encode('ascii')
    file2 = fopen(pipe2_bytes, 'r')

    output_read_file_bytes = output_read_file.encode('ascii')
    out = fopen(output_read_file_bytes, 'w')
    
    umi_indices = dict([seq, ind] for ind, seq in
                       enumerate(util.get_all_kmers(6)))    
    
    c = 0
    w1_len = len(w1)
    passed_w1 = 0
    passed_bc = 0
    passed_qual = 0
    t0 = time.time()
    
    while True:
        if c_max_reads !=0 and c >= max_reads:
            break

        # read barcode read name (which we ignore)
        ret = fgets(buf, buf_size, file1)
        if ret == NULL:
            # end of file reached, so we're done!
            break
            
        c+=1

        # read barcode read sequence
        fgets(seq, buf_size, file1)
        r1 = seq.decode('ascii')
        r1_len = strlen(seq)
        
        # read remaining 2 lines of barcode read block
        fgets(buf, buf_size, file1)    
        fgets(qual, buf_size, file1)

        # read the four lines of mRNA read block
        fgets(buf, buf_size, file2)
        fgets(l2, buf_size, file2)
        fgets(l3, buf_size, file2)
        fgets(l4, buf_size, file2)

        # try to identify W1 sequence
        bc1_len = 0
        if r1_len >= 48:  # make sure read is long enough
            if r1[8:(8+w1_len)] in w1_mapping:
                bc1_len = 8
            elif r1[9:(9+w1_len)] in w1_mapping:
                bc1_len = 9
            elif r1[10:(10+w1_len)] in w1_mapping:
                bc1_len = 10
            elif r1[11:(11+w1_len)] in w1_mapping:
                bc1_len = 11

        if bc1_len == 0:
            # we couldn't identify the W1 sequence, so we'll skip this read
            continue
            
        passed_w1 += 1

        bc1_found = 1
        bc2_found = 1

        # try to identify first barcode
        try:
            bc1 = bc1_exact_mapping[r1[:bc1_len]]
        except KeyError:
            try:
                bc1 = bc1_mismatch_mapping[r1[:bc1_len]]
            except KeyError:
                bc1_found = 0

        bc2_start = bc1_len + w1_len

        # try to identify second barcode
        try:
            bc2 = bc2_exact_mapping[r1[bc2_start:(bc2_start+8)]]
        except KeyError:
            try:
                bc2 = bc2_mismatch_mapping[r1[bc2_start:(bc2_start+8)]]
            except KeyError:
                bc2_found = 0

        if bc1_found == 0 or bc2_found == 0:
            # we failed to identify at least one barcode, so we'll skip this
            # read
            continue
                
        passed_bc += 1
                
        # extract UMI sequence
        umi_start = bc2_start + 8
        umi = r1[umi_start:(umi_start+6)]

        # check whether UMI sequence contains an "N"
        good_qual = 1
        for p in range(umi_start, umi_start+6):
            if seq[p] == 78:
                good_qual = 0
        
        if good_qual == 0:
            # UMI sequence contained an "N", so we'll skip this read
            continue

        passed_qual += 1

        bc1_counts[bc1] += 1
        bc2_counts[bc2] += 1
        bc_counts[bc1, bc2] += 1


        # look up the UMI index
        umi_idx = umi_indices[umi]

        # extract the UMI quality scores
        for p in range(umi_start, umi_start+6):
            # we use PHRED+66 to avoid problems with read names (e.g., "/")
            qual[p] = qual[p] + 33
        qual_str = qual.decode('ascii')[umi_start:(umi_start+6)]

        # create the read name
        r1 = '@%03d-%03d_%04d_%s_%d\n' \
                % (bc1, bc2, umi_idx, qual_str, c)

        # write the read to the output FASTQ file
        r1_bytes = r1.encode('ascii')  # needs to be assigned a Python variable
        fputs(<char*>r1_bytes, out)
        fputs(l2, out)
        fputs(l3, out)
        fputs(l4, out)                

    fclose(out)
    fclose(file1)
    fclose(file2)
    os.remove(pipe1)
    os.remove(pipe2)
    
    t1 = time.time()
    t = t1-t0
    
    _LOGGER.info('Processed %d reads in %.1f s.', c, t)
    _LOGGER.info("(That's %.1f s per million reads.)", 1e6 *(t/c))
    _LOGGER.info('Reads with W1 sequence found: %d / %d (%.1f %%)',
                 passed_w1, c, 100*(passed_w1/float(c)))
    _LOGGER.info('- of which with valid barcodes: %d / %d (%.1f %%)',
                 passed_bc, passed_w1, 100*(passed_bc/float(passed_w1)))
    _LOGGER.info('-- of which with valid UMI sequences: %d / %d (%.1f %%)',
                 passed_qual, passed_bc, 100*(passed_qual/float(passed_bc)))        
    _LOGGER.info('Total fraction of reads passed: %d / %d (%.1f %%)',
                 passed_qual, c, 100*(passed_qual/float(c)))

    bc1_index = pd.Index(barcodes1)
    bc2_index = pd.Index(barcodes2)
    df = pd.DataFrame(np.uint64(bc_counts), index=bc1_index, columns=bc2_index)
    df.to_csv(output_count_file, sep='\t')    
    
    _LOGGER.info('Done.')
