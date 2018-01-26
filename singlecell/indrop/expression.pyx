#cython: language_level=3, boundscheck=False
"""Functions for quantifying gene expression based on inDrop data."""

import logging
import time
import os
import csv
from collections import Counter

import pandas as pd
from scipy import sparse
import numpy as np
cimport numpy as np
np.import_array()

from genometools.expression import ExpMatrix, ExpGeneTable

from .. import util

_LOGGER = logging.getLogger(__name__)


def quantify_expression(
        barcode_counts_file,
        chromosome_length_file,
        read_info_dir,
        gene_file, genome_annotation_file,
        num_cells,
        gene_expression_output_file, transcript_expression_output_file,
        min_umi_qual=20, cell_prefix='cell_',
        dense_gene_expression_output_file=None,
        max_chroms=None):
    """Determines UMIFM reads per gene for inDrop data.
    
    TODO: docstring"""
    
    cdef int n, p, j, k
    cdef unsigned short max_gene_index, gene_idx
    cdef long long int cur, cur2, pos, num_transcripts_seen    
    cdef long long int[::1]  bc_gene
    cdef unsigned short[::1] all_genes
    cdef unsigned short[::1] gene_max_transcripts_view
    cdef unsigned short[:, ::1] cell_indices
    cdef int bc1_idx, bc2_idx
    
    if num_cells >= 65536:
        raise ValueError('The max. number of cells supported is 65535!')

    min_umi_qual_ascii = min_umi_qual + 66
                
    # read chromosome lengths
    chromlen = pd.read_csv(chromosome_length_file,
                           sep='\t', index_col=0, squeeze=True)
    chromlen_dict = dict([c, l] for c, l in chromlen.iteritems())

    # read gene list
    gene_table = ExpGeneTable.read_tsv(gene_file)

    # get gene names that are guaranteed to be unique
    gene_names = util.get_readable_gene_identifiers(gene_table)

    # series with index = Ensembl ID, value = unique gene name
    genes = pd.Series(index=gene_table.index, data=gene_names)
    
    if len(gene_names) > 65536:
        raise ValueError('The max. number of genes supported is 65536!')
    
    gene_indices = dict([name, i] for i, name in enumerate(gene_names))

    # create dictionary mapping chromosomes to a list of genes on those chrom.
    chrom_genes = dict([chrom, []] for chrom in chromlen.index)
    for ensembl_id, name in genes.items():
        chrom_genes[gene_table.loc[ensembl_id, 'chromosome']].append(name)

    # create dictionary mapping chromosomes to numpy ndarray of gene indices
    chrom_gene_indices = dict()
    for chrom in chromlen.index:
        chrom_gene_indices[chrom] = \
                np.uint16([gene_indices[name]
                           for name in chrom_genes[chrom]])
    
    # read barcode counts and sort
    counts = pd.read_csv(barcode_counts_file, sep='\t', header=0, index_col=0)
    a = np.argsort(counts.values.ravel())[::-1]
    x, y = np.unravel_index(a, counts.shape)
    
    # generate matrix that contains sample numbers
    # (defining which barcode combinations we'll include)
    cell_indices = 65535 * np.ones((384, 384), dtype=np.uint16)
    cell_names = ['']*num_cells
    for i in range(num_cells):
        cell_indices[x[i], y[i]] = i
        cell_names[i] = '%s%03d-%03d' % (cell_prefix, x[i], y[i])

    p = len(gene_names)
    n = num_cells

    _LOGGER.info('Total number of protein-coding genes: %d', p)
    
    num_total = counts.values.ravel().sum()
    num_included = counts.values.ravel()[a[:n]].sum()
    
    # intialize "gene expression" array
    # (holds the number of unique umis per gene in each cell)
    G = np.zeros((p, n), dtype=np.uint16)
    
    transcript_arrays = []
    transcript_names = []
    
    num_total = 0
    num_good_qual = 0
    num_included = 0
    gene_max_transcripts = np.zeros(p, dtype=np.uint16)
    gene_max_transcripts_view = gene_max_transcripts
    t0 = time.time()
    for ic, chrom in enumerate(chromlen.index):

        if max_chroms is not None and ic >= max_chroms:
            break

        read_info_file = os.path.join(read_info_dir, '%s.tsv' % chrom)        
        if not os.path.isfile(read_info_file):
            _LOGGER.warning('Read info file for chromosome "%s" not found. '
                            'Skipping...', chrom)
            continue

        # get number of reads from file header
        l = open(read_info_file).readline()
        assert l.startswith('#num_reads=')
        num_reads = int(l[11:-1])

        # => each row in the array: (umi index, gene index, sample index)
        read_info = np.empty((num_reads, 3), dtype=np.uint16)

        i = 0
        with open(read_info_file, encoding='ascii', newline='') as fh:
            reader = csv.reader(fh, dialect='excel-tab')
            next(reader)  # skip header
            for l in reader:
                num_total += 1
                
                # figure out which cell the read came from
                barcode = l[0]
                bc1_idx = int(barcode[:3])  # first barcode index
                bc2_idx = int(barcode[4:7])  # second barcode index
                cell_idx = cell_indices[bc1_idx, bc2_idx]

                if cell_idx < n:
                    num_included += 1

                    # check if UMI sequence quality is OK
                    if min(l[2].encode('ascii')) >= min_umi_qual_ascii:
                        num_good_qual += 1
                        
                        target_gene = l[3]
                        gene_idx = gene_indices[target_gene]
                        umi_idx = int(l[1])
                        read_info[i, 0] = umi_idx
                        read_info[i, 1] = gene_idx
                        read_info[i, 2] = cell_idx
                        i += 1

        # we probably excluded some reads -
        # we selected specific cellular barcodes,
        # and we also filtered for UMI sequence quality -
        # so we have to truncate the array
        read_info = read_info[:i, :]

        # we now sort the array using cell index as primary key,
        # gene index as secondary key,
        # and UMI index as tertiary key
        a = np.lexsort(read_info.T)
        read_info = read_info[a, :]

        # count the number of reads per cell
        bc_cell = np.bincount(read_info[:, 2])

        ### quantify gene expression

        # navigate the sorted read array to count UMIs
        # for each gene in each cell
        cur = 0
        for j in np.nonzero(bc_cell)[0]:
            # count the number of reads per gene for this cell
            bc_gene = np.bincount(read_info[cur:(cur+bc_cell[j]), 1])
            cur2 = 0
            for k in np.nonzero(bc_gene)[0]:
                # determine number of unique UMIs
                # for all reads from this gene (for the current cell)
                G[k, j] = np.unique(
                    read_info[(cur+cur2):(cur+cur2+bc_gene[k]), 0]).size
                # go to the reads from the next gene
                cur2 += bc_gene[k]
            # go to the reads from the next cell
            cur += bc_cell[j]

        ### quantify transcript expression

        # determine the max. number of transcripts across all cells,
        # for each gene on this chromosomes
        sel = chrom_gene_indices[chrom]
        gene_max_transcripts[sel] = G[sel, :].max(axis=1)
        #genes_seen = chrom_gene_indices[chrom_gene_max_transcripts > 0]

        # the sum of all those values determines the number of rows
        # required for the matrix holding the number of reads per
        # expressed transcript (from all genes on this chromosome)
        num_transcript_rows = gene_max_transcripts[sel].sum()
        T = np.zeros((num_transcript_rows, n), dtype=np.uint32)

        # determine all genes (on this chromosome)
        # with at least one transcript
        # note: np.unique returns a sorted array
        all_genes = np.unique(read_info[:, 1])

        # navigate the sorted read array to count UMIs for each
        # transcript in each cell
        cur = 0
        for j in np.nonzero(bc_cell)[0]:

            # count the number of reads per gene for this cell
            bc_gene = np.bincount(read_info[cur:(cur+bc_cell[j]), 1],
                                  minlength=p)

            # go over ALL genes expressed on this chromosome
            # => we have to go over genes in the same order for each cell,
            # in order to be able to store expression values in the right
            # rows
            cur2 = 0
            pos = 0
            for k in range(all_genes.size):
                gene_idx = all_genes[k]
                # if this gene is expressed in this cell,
                # we quantify the expression of its transcripts
                if bc_gene[gene_idx] > 0:
                    # count reads for each UMI present
                    counts = Counter(
                        read_info[(cur+cur2):
                                  (cur+cur2+bc_gene[gene_idx]), 0])
                    # sort present UMIs by their number of reads
                    sorted_umis = sorted(counts.keys(),
                            key=lambda x: -counts[x])
                    T[pos:(pos+len(counts)), j] = \
                            [counts[u] for u in sorted_umis]
                    # go to the reads from the next gene
                    cur2 += bc_gene[gene_idx]

                # move down the approriate number of rows
                # in the transcript expression array (T)
                pos += (<unsigned long long int>
                        (gene_max_transcripts_view[gene_idx]))
            # go to the reads from the next cell
            cur += bc_cell[j]

        chrom_transcript_names = \
                ['%s_%d' % (gene_names[i], j)
                 for i in all_genes
                 for j in range(gene_max_transcripts_view[i])]

        transcript_names.append(chrom_transcript_names)
        transcript_arrays.append(T)
    
    T = np.concatenate(transcript_arrays, axis=0)
    transcript_names = [item for subl in transcript_names for item in subl]
    assert len(transcript_names) == T.shape[0]
    
    t1 = time.time()
    t = t1-t0
    _LOGGER.info('Processed %d reads overlapping exon(s) from a single gene in %.1f s.', num_total, t)
    _LOGGER.info("(That's %.1f s per million aligned reads.)",
                 1e6*(t/float(num_total)))            
    
    _LOGGER.info('Number of reads from included cells: %d (%.1f%%)',
                 num_included, 100*(num_included/float(num_total)))
    if num_included > 0:
        _LOGGER.info('- of which with good UMI sequence quality: %d (%.1f%%)',
                     num_good_qual, 100*(num_good_qual/float(num_included)))
    
    # write sparse gene expression output file
    matrix = ExpMatrix(X=G, genes=gene_names, cells=cell_names)
    matrix.sort_index(inplace=True)
    matrix.write_sparse(gene_expression_output_file)
    
    if dense_gene_expression_output_file and \
            dense_gene_expression_output_file != gene_expression_output_file:
        # write dense gene expression output file
        matrix.write_tsv(dense_gene_expression_output_file)

    ### write sparse transcript expression output file
    matrix = ExpMatrix(X=T, genes=transcript_names, cells=matrix.cells)
    matrix.sort_index(inplace=True)
    matrix.write_sparse(transcript_expression_output_file)
