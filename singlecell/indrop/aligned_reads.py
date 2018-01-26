"""Functions for processing aligned reads."""

import logging
from math import ceil
import time
import multiprocessing
import os
import tempfile
import shutil
from collections import Counter, OrderedDict

import pandas as pd
import pysam
import HTSeq

from genometools import misc
from genometools.expression import ExpGeneTable
from .. import util

_LOGGER = logging.getLogger(__name__)
_ga = None  # global HTSeq.GenomicArrayOfSets


def _get_chromosome_chunks(length, num_jobs):
    """Returns evenly sized intervals (start, end) covering the chromosome.
    
    Parameters
    ----------
    length : int
        The length of the chromosome.
    num_jobs : int
        The number of intervals.

    Returns
    -------
    list of (start, end) tuples
      The interval start and end positions.      
    """
    chunk_size = int(ceil((length/float(num_jobs))))
    chunks = []
    for i in range(num_jobs):
        chunks.append((i*chunk_size, min((i+1)*chunk_size, length)))
    return chunks


def _process_chromosome(alignment_file, chromosome, start, end, output_file):
    """Process the aligned reads covering a part (or all) of a chromosome.

    The output is a tab-delimited text file where each row corresponds to one
    read, and the columns contain its cellular barcode, its UMI barcode, its
    UIM barcode sequencing qualities, and the gene it was mapped to,
    respectively.

    Parameters
    ----------
    alignment_file : str
        Path of the alignment (.bam) file. Must be accompanied by an index
        (.bai) file.
    chromosome : str
        The name of the chromosome.
    start : int
        The start position of the interval to be processed (0-indexed).
    end : int
        The end position of the interval to be processed (exclusive).
    output_file : str
        Path of the output (.tsv) file.
    
    Returns
    -------
    str
        The chromosome name.
    int
        The total number of reads covering the interval.
    int
        The number of reads overlapping one or more exons from a single gene.
    int
        The number of reads overlapping exons from multiple genes.

    Notes
    -----
    This function relies on global variable storing a
    `HTSeq.GenomicArrayOfSets` instance. This is to reduce the overhead when
    submitting a call of this function to a `multiprocessing.Pool` instance.
    """
    num_reads = 0
    num_unique = 0
    num_ambig = 0
    with pysam.AlignmentFile(alignment_file, 'rb') as bam, \
            open(output_file, 'w', encoding='ascii') as ofh:
        for read in bam.fetch(chromosome, start, end):
            num_reads += 1
            read_name = read.query_name

            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            iv_seq = (HTSeq.GenomicInterval(chromosome, b[0], b[1], strand)
                      for b in read.get_blocks())
            fs = set()
            for iv in iv_seq:
                for iv2, fs2 in _ga[iv].steps():
                    fs = fs.union(fs2)

            if len(fs) > 1:
                # the aligned read overlaps with exons from more than one gene
                # => we exclude those reads
                num_ambig += 1

            elif len(fs) == 1:
                # the aligned read overlaps with one or more exons
                # from the same gene
                num_unique += 1
                barcode = read_name[:7]
                umi_idx = read_name[8:12]
                umi_qual = read_name[13:19]
                target_gene = list(fs)[0]                
                ofh.write('%s\t%s\t%s\t%s\n' % (barcode, umi_idx, umi_qual, target_gene))
                
    return chromosome, num_reads, num_unique, num_ambig


def process_aligned_reads(
        alignment_file, chromosome_length_file,
        gene_file, gene_name_file, genome_annotation_file, output_dir,
        num_jobs=1):
    """Process aligned inDrop sequencing reads.
    
    
    """
    
    global _ga

    t0 = time.time()

    # make sure BAM file is readable and has an index
    with pysam.AlignmentFile(alignment_file, 'rb') as bam:
        if not bam.check_index():
            raise RuntimeError('The bam file is missing an index.')

    # make sure output directory exists
    misc.make_sure_dir_exists(output_dir)

    # read chromosome lengths
    chromlen = pd.read_csv(chromosome_length_file,
                           sep='\t', index_col=0, squeeze=True)

    # read gene list
    gene_table = ExpGeneTable.read_tsv(gene_file)

    # intialize HTSeq `GenomicArrayOfSets` object for storing exons
    chromlen_dict = dict([c, l] for c, l in chromlen.iteritems())   
    _ga = HTSeq.GenomicArrayOfSets(chromlen_dict, stranded=True)
    
    # parse GTF file and extract all exons
    #logger = logging.getLogger('singlecell.util')
    #logger.setLevel(logging.ERROR)  # disable logging from get_gene_exons
    gene_exons = util.get_gene_exons(gene_table, genome_annotation_file)
    #logger.setLevel(logging.INFO)

    # get gene names that are guaranteed to be unique
    #gene_names = pd.Series(
    #    index=gene_table.index,
    #    data=util.get_readable_gene_identifiers(gene_table))
    gene_names = pd.read_csv(gene_name_file, sep='\t',
                             index_col=0, squeeze=True)

    # create HTSeq.GenomicArrayOfSets containing the exons
    for id_, name in gene_names.items():
        gene = gene_table.loc[id_]
        chrom = gene['chromosome']
        exons = gene_exons[id_]
        merged = util.merge_intervals(exons)

        if gene['position'] >= 0:  # "+" strand orientation
            strand = '+'
        else:
            strand = '-'

        for m in merged:
            iv = HTSeq.GenomicInterval(chrom, m[0], m[1], strand)
            _ga[iv] += name
    
    t1 = time.time()
    t = t1 - t0
    _LOGGER.info('Preparation time: %.2f s.', t)
    
    # sort chromosomes by their length
    # (we want the mitochondrial chromosome to be among
    #  the first submitted, since it can potentially have
    #  a lot of reads, but it doesn't lend itself to
    #  chunking very well - similarly for ribosomal
    #  contigs)
    chromosomes = sorted(chromlen_dict.keys(),
                            key=lambda c:chromlen_dict[c])

    temp_dir = tempfile.mkdtemp(prefix='indrop_')
    _LOGGER.debug('Temp dir: %s', temp_dir)
            
    pool = multiprocessing.Pool(num_jobs)
    
    chrom_chunks = dict()
    chrom_chunk_files = dict()
    chrom_reads = dict()
    jobs = []
    t0 = time.time()
    for j, chrom in enumerate(chromosomes):
        
        # only split up chromosome into chunks if it is large (>= 1M bp)
        clen = chromlen_dict[chrom]            
        if num_jobs == 1 or clen < 1000000:     
            chunks = [(0, clen)]
        else:
            chunks = _get_chromosome_chunks(clen, num_jobs)
            assert len(chunks) == num_jobs

        chunk_files = [
            os.path.join(temp_dir, '%s_%02d.tsv' % (chrom, i))
            for i in range(len(chunks))]
        chrom_chunk_files[chrom] = chunk_files
        
        # remember the number of chunks for this chromosome
        chrom_chunks[chrom] = len(chunks)

        # submit processing job(s) to the pool
        for i, ch in enumerate(chunks):
            jobs.append(pool.apply_async(
                _process_chromosome,
                (alignment_file, chrom, ch[0], ch[1],
                    chunk_files[i])))

    t1 = time.time()
    t = t1 - t0
    _LOGGER.info('Submitted a total of %d jobs in %.1f s.', len(jobs), t)

    # wait for all threads to finish
    t0 = time.time()
    pool.close()  # we need to call this before pool.join()
    pool.join()  # waits for all jobs to be done
    t1 = time.time()

    # count up reads (total/unique/ambiguous)
    num_reads = 0
    num_unique = 0
    num_ambig = 0
    chrom_reads = Counter()
    for j in jobs:
        # make sure job is done and had no error
        assert j.ready() and j.successful()
        res = j.get()  # get the function return values
        chrom = res[0]
        num_reads += res[1]
        num_unique += res[2]
        num_ambig += res[3]
        # also count number of unique reads for each chromosome separately
        chrom_reads[chrom] += res[2]  

    # combine chunked chromosome files
    for chrom in chromosomes:
        output_file = os.path.join(output_dir, '%s.tsv' % chrom)
        # write a header line containing the number of reads in this file
        with open(output_file, 'w') as ofh:
            ofh.write('#num_reads=%d\n' % chrom_reads[chrom])
        util.concatenate_files(
            chrom_chunk_files[chrom], output_file, append=True)

    t = t1-t0
    _LOGGER.info('Processed %d reads in %.1f s. That\'s %.2f s per million '
                 'reads.',
                 num_reads, t, 1e6*(t/float(num_reads)))

    _LOGGER.info('The number of unique reads is: %d (%.1f %%)',
                 num_unique, 100*(num_unique/float(num_reads)))
    _LOGGER.info('The number of ambiguous reads is: %d (%.1f %%)',
                 num_ambig, 100*(num_ambig/float(num_reads)))
    

    shutil.rmtree(temp_dir)
