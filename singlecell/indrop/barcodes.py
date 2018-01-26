"""Functions for working with and analyzing inDrop barcodes."""

from collections import Counter, defaultdict
import logging
import time
import sys

import pandas as pd
import numpy as np
import pysam

from plotly.offline import plot
import plotly.graph_objs as go

from genometools.expression import ExpMatrix, ExpGeneTable

from .. import util
from .barcodes_cython import count_mapped_reads

_LOGGER = logging.getLogger(__name__)


def get_exact_mapping(barcode_file):
    """
    Returns a mapping from sequences to barcode indices, for
    all sequences matching barcodes exactly.
    """
    # read the barcodes
    barcodes = pd.read_csv(barcode_file, squeeze=True, header=None)

    # take reverse complement sequence
    barcodes = barcodes.apply(util.get_reverse_complement)

    mapping = dict([bc, i] for i, bc in enumerate(barcodes))
    return mapping


def get_mismatch_mapping(barcode_file):
    """
    Returns a mapping from sequences to barcode indices, for
    all sequences with one mismatch, excluding abiguous cases.
    """
    # read the barcodes
    barcodes = pd.read_csv(barcode_file, squeeze=True, header=None)

    # take reverse complement sequence
    barcodes = barcodes.apply(util.get_reverse_complement)

    # exact matches should always be preferred to 1-mismatch matches
    exact = set(barcodes)
    assert len(exact) == len(barcodes)

    # generate a list of all 1-mismatch mathces
    # to do so, generate a dictionary of "mismatch sequence"
    #   => ["barcode sequence 1", ...]
    # this lets us map back to how many mismatches were generated for each
    # barcode

    mismatch = defaultdict(set)
    for i, bc in enumerate(barcodes):
        for mm in util.get_mismatch_sequences(bc):
            if mm not in exact:
                mismatch[mm].add((bc, len(bc), i))

    _LOGGER.debug('How many barcodes are assigned to each mismatch sequence?')
    _LOGGER.debug(str(Counter(len(v) for v in mismatch.values())))

    bc_count = Counter()
    for v in mismatch.values():
        for bc in v:
            bc_count[bc[0]] += 1
    _LOGGER.debug('In how many mismatch sequences does each barcode occur?')
    _LOGGER.debug(str(Counter(bc_count.values())))

    # only keep mismatch mapping
    mismatch_unambiguous = dict(
        (mm, list(barcodes)[0][2])
        for mm, barcodes in mismatch.items() if len(barcodes) == 1)
    _LOGGER.debug('Found %d unambiguously mapping mismatch sequences.',
                  len(mismatch_unambiguous))

    return mismatch_unambiguous


def count_transcriptomic_reads(
            alignment_file, barcode1_file, barcode2_file,
            gene_file, genome_annotation_file,
            output_file, max_genes=None):
    """Count the number of transcriptomic reads for each inDrop barcode.
    
    TODO: docstring"""
    
    counts = np.zeros((384, 384), dtype=np.uint64)

    genome = ExpGeneTable.read_tsv(gene_file)
    p = genome.num_genes
    _LOGGER.info('Number of protein-coding genes: %d', p)
    
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
    
    # only parse GTF after we've made sure BAM file is readable and has an
    # index
    with pysam.AlignmentFile(alignment_file, 'rb') as bam:
        if not bam.check_index():
            raise AssertionError('Bam file is missing the index!')

    gene_exons = util.get_gene_exons(genome, genome_annotation_file)    
    gene_names = sorted(list(gene_exons.keys()))
    
    total = 0
    valid = 0
    t0 = time.time()
    with pysam.AlignmentFile(alignment_file, 'rb') as bam:

        for i, gene_name in enumerate(gene_names):
  
            if (i+1) % 1000 == 0:
                print(i+1, end=' ')
                sys.stdout.flush()
                
            if max_genes is not None and i+1 > max_genes:
                break
  
            exons = gene_exons[gene_name]
            merged = util.merge_intervals(exons)
            
            gene = genome[gene_name]
            if gene.position >= 0:
                filter_func = lambda x: not x.is_reverse
            else:
                filter_func = lambda x: x.is_reverse            
            
            for iv in merged:
                for read in bam.fetch(gene.chromosome, iv[0], iv[1]):
                    total += 1
                    if filter_func(read):
                        valid += 1
                        bc1 = int(read.query_name[:3])  # first barcode index
                        bc2 = int(read.query_name[4:7])  # second barcode index                        
                        counts[bc1, bc2] += 1
    t1 = time.time()
                    
    print()
    t = t1-t0
    _LOGGER.info('Processed %d alignments in %.1f s.', total, t)
    _LOGGER.info("(That's %.1f s per million reads.)", 1e6*(t/total))            
    _LOGGER.info('Reads with valid orientation: %d (%.1f%%)'
                 % (valid, 100*(valid/float(total))))
    
    bc1_index = pd.Index(barcodes1)
    bc2_index = pd.Index(barcodes2)
    df = pd.DataFrame(counts, index=bc1_index, columns=bc2_index)
    df.to_csv(output_file, sep='\t')


def plot_read_count_distribution(barcode_count_file, output_file,
                                 xaxis_label=('# mapped reads '
                                              '(log<sub>10</sub>-scale)')):
    """Plot histogram of the distribution of reads per barcode.
    
    TODO: docstring""" 
    matrix = ExpMatrix.read_tsv(barcode_count_file)

    x = np.float64(matrix.values.ravel())
    num_total_reads = int(np.sum(x))
    x[x < 1] = 1
    x = np.log10(x)

    data = [
        go.Histogram(x=x, nbinsx=100)
    ]

    layout = go.Layout(
        title='Total number of mapped reads: %d' % num_total_reads,
        font=dict(
            size=20,
            family='serif',
        ),
        xaxis=dict(
            title=xaxis_label,
        ),
        yaxis=dict(
            title='# barcodes',
            type='log',
        ),
    )

    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename=output_file, show_link=False, auto_open=False)
