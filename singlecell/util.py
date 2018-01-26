"""Utility functions."""

import subprocess
import logging
import os
import shutil
import stat
import itertools
from collections import OrderedDict
from pkg_resources import resource_string

import pandas as pd
from genometools.expression import ExpGeneTable
from genometools import gtf

import singlecell


_LOGGER = logging.getLogger(__name__)


def get_readable_gene_identifiers(gene_table: ExpGeneTable):
    """Return unique gene identifiers that primarily use the genes' names."""
    # count occurrences for each of gene name
    counts = gene_table['name'].value_counts()
    gene_counts = counts.loc[gene_table['name']]
    gene_ids = gene_table.index.tolist()
    
    gene_ids = [name if c == 1 else '%s_%s' % (name, gene_ids[i])
                for i, (name, c) in enumerate(gene_counts.items())]

    return gene_ids


def get_edit_sequences(seq, num_edits, bases=None):
    """Return all nucleotide sequences with a given hamming distance."""

    if num_edits > len(seq):
        raise ValueError('Asked to make make more edits (%d) than the length '
                         'of the sequence (%d nt).' % (num_edits, len(seq)))

    if bases is None:
        bases = set('ACGT')

    length = len(seq)
    all_bases = [bases for i in range(num_edits)]
    seq_list = [nt for nt in seq]
    mismatch = []   
    for comb in itertools.combinations(range(length), num_edits):
        for subs in itertools.product(*all_bases):
            mut = seq_list[:]
            valid = True
            for pos, nt in zip(comb, subs):
                if mut[pos] == nt:
                    valid = False
                    break
                mut[pos] = nt
            if valid:
                mismatch.append(''.join(mut))
    
    return sorted(mismatch)


def concatenate_files(input_files, output_file, append=False):
    write_mode = 'wb'
    if append:
        write_mode = 'ab'
    with open(output_file, write_mode) as ofh:
        for f in input_files:
            with open(f, 'rb') as ifh:
                shutil.copyfileobj(ifh, ofh, 16*1024*1024)


def make_file_executable(path):
    """Sets the user executable flag for a file."""
    st = os.stat(path)
    os.chmod(path, st.st_mode | stat.S_IEXEC)


def zcat_subproc(path):
    """Creates a subprocess for decompressing a gzip file.
    
    TODO: docstring"""
    subproc = subprocess.Popen('gunzip -c "%s"' % path, shell=True,
                               stdout=subprocess.PIPE)
    return subproc


def get_all_kmers(k, kmer='', kmer_list=None):
    """Returns all possible k-mer sequences (for A/C/G/T alphabet).
    
    TODO: docstring"""
    
    if kmer_list is None:
        kmer_list = []
    if len(kmer) == k:
        kmer_list.append(kmer)
    else:
        for nuc in ['A', 'C', 'G', 'T']:
            var = kmer + nuc
            get_all_kmers(k, var, kmer_list)
    if not kmer:
        return kmer_list


def get_mismatch_sequences(seq):
    """Generates all nucleotide sequences with hamming distance 1 to `seq`.
    
    TODO: docstring"""

    for pos in range(len(seq)):
        for nuc in ['A', 'C', 'G', 'T']:
            if nuc != seq[pos]:
                mm = seq[:pos] + nuc + seq[(pos+1):]
                yield mm


def get_reverse_complement(seq):
    """Returns the reverse complement of a nucleotide sequence.
    
    TODO: docstring"""
    rc = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    compseq = ''.join([rc[nuc] for nuc in seq[::-1]])
    return compseq


def get_gene_exons(gene_table, genome_annotation_file, chunksize=10000):
    """Parse GTF file and get a dictionary of gene=>list of exon intervals.
    
    (Only for protein-coding genes.)
    TODO: docstring"""
    
    # get gene names that are guaranteed to be unique
    #gene_names = get_readable_gene_identifiers(gene_table)

    # series with index = Ensembl ID, value = unique gene name
    #genes = pd.Series(index=gene_table.index, data=gene_names)

    # sort genes by chromosome, strand, and then position
    sorted_gene_ids = sorted(
        [id_ for id_ in gene_table.index],
        key=lambda id_: [gene_table.loc[id_, 'chromosome'],
                         gene_table.loc[id_, 'position'] < 0,
                         abs(gene_table.loc[id_, 'position'])])
    #genes = genes.loc[sorted_gene_ids]
    gene_table = gene_table.loc[sorted_gene_ids]

    # dictionary for holding list of intervals for each gene
    gene_exons = OrderedDict([id_, []] for id_ in gene_table.index)

    valid = 0
    total = 0
    
    _LOGGER.info('Parsing GTF file "%s" in chunks...', genome_annotation_file)

    for i, df in enumerate(pd.read_csv(
            genome_annotation_file, dtype={0: str},
            sep='\t', comment='#', header=None, chunksize=chunksize)):

        # select only exon entries
        df_sel = df.loc[df.iloc[:, 2] == 'exon']

        # extract gene IDs
        gene_ids = df_sel.iloc[:, 8].apply(
            lambda x: gtf.parse_attributes(x)['gene_id'])

        for id_, chrom, start, end in zip(
                gene_ids,
                df_sel.iloc[:, 0], df_sel.iloc[:, 3], df_sel.iloc[:, 4]):

            total += 1
            
            try:
                gene = gene_table.loc[id_]
            except KeyError:
                # this gene is not contained in the gene table
                continue

            gene_chrom = gene_table.loc[id_, 'chromosome']
            if chrom != gene_chrom:
                _LOGGER.warning('%s exon ignored (wrong chromosome: '
                                '%s instead of %s).',
                                id_, chrom, gene_chrom)
            else:
                valid += 1
                gene_exons[id_].append([start-1, end])

    _LOGGER.info('%d / %d exons from valid genes (%.1f %%).',
                 valid, total, 100*(valid/float(total)))
    
    return gene_exons


def merge_intervals(intervals):
    """Merge overlapping intervals.
    
    TODO: docstring"""
    
    if not intervals:
        return []
    # sort intervals by start position
    intervals = sorted(intervals, key=lambda x:x[0])
    merged = []

    cur = list(intervals[0])
    for iv in intervals[1:]:
        # interval starts inside/right after current interval
        if iv[0] <= cur[1]:
            if iv[1] > cur[1]:  # interval ends after current interval
                cur[1] = iv[1]
        else:
            merged.append(cur)
            cur = list(iv)
    merged.append(cur)
    return merged


def get_mitochondrial_genes(species='human'):
    """Get a list of all mitochondrial genes for a given species.
    
    "Mitochondrial genes" are defined here as all genes on the mitochondrial
    chromosome.

    TODO: docstring
    """
    path = os.path.join(singlecell._root,
                        'data', 'gene_lists', 'mitochondrial_%s.tsv' % species)
    with open(path) as fh:
        return fh.read().split('\n')


def get_ribosomal_genes(species='human'):
    """Get a list of all ribosomal genes for a given species.
    
    "Ribosomal genes" are defined here as all protein-coding genes whose
    protein products are a structural component of the small or large ribosomal
    subunit (including fusion genes).

    TODO: docstring
    """

    path = os.path.join(singlecell._root,
                        'data', 'gene_lists', 'ribosomal_%s.tsv' % species)
    with open(path) as fh:
        return fh.read().split('\n')


def get_plotly_js():
    """Return the plotly javascript code.
    
    TODO: docstring
    """
    # resource_string?
    path = 'package_data/plotly.min.js'
    return resource_string('plotly', path).decode('utf-8')


def is_empty_dir(dir_):
    """Tests whether a directory is empty.

    Note: Also returns True if the directory doesn't exist.

    TODO: docstring
    """
    is_empty = True
    try:
        _, dirnames, filenames = next(os.walk(dir_))
        if dirnames or filenames:
            is_empty = False
    except StopIteration:
        pass
    return is_empty
