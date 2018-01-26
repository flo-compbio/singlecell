"""Fixtures for the inDrop tests."""

import logging
import pathlib

import pytest

from genometools import misc

_LOGGER = logging.getLogger(__name__)


@pytest.fixture(scope='session')
def my_output_pypath(tmpdir_factory):
    """temporary directory for storing test data"""
    pypath = tmpdir_factory.mktemp('singlecell_output')
    return pypath


@pytest.fixture(scope='session')
def my_barcode_read_file(my_data_pypath):
    """barcode reads"""
    _LOGGER.info('Downloading the mRNA read file...')
    url = 'https://www.dropbox.com/s/has2sys148jwlwf/SRR3879606_10k_1.fastq.gz?dl=1'
    path = str(my_data_pypath.join('barocde_reads.fastq.gz'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_mrna_read_file(my_data_pypath):
    """mRNA reads"""
    _LOGGER.info('Downloading the mRNA read file...')
    url = 'https://www.dropbox.com/s/nxu3o28mj7kpm4l/SRR3879606_10k_2.fastq.gz?dl=1'
    path = str(my_data_pypath.join('mRNA_reads.fastq.gz'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_alignment_file(my_data_pypath):
    """Alignment file."""
    _LOGGER.info('Downloading the alignment file (+ index file)...')
    
    url = 'https://www.dropbox.com/s/1bkwa3puy504eys/SRR3879606_10k.bam.bai?dl=1'
    path = str(my_data_pypath.join('alignment.bam.bai'))
    misc.http_download(url, path)
    
    url = 'https://www.dropbox.com/s/jx0049j7w6pbei3/SRR3879606_10k.bam?dl=1'
    path = str(my_data_pypath.join('alignment.bam'))
    misc.http_download(url, path)

    return path

@pytest.fixture(scope='session')
def my_chromosome_length_file(my_data_pypath):
    """Chromosome length file."""
    url = 'https://www.dropbox.com/s/mw8o9ur48pr1squ/chromosome_lengths_human_GRCh38.tsv?dl=1'
    path = str(my_data_pypath.join('chromosome_lengths_human_GRCh38.tsv'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_gene_file(my_data_pypath):
    """Protein-coding gene file."""
    _LOGGER.info('Downloading the gene file...')
    url = 'https://www.dropbox.com/s/o3tfaubrcoe92cr/protein_coding_genes_human_ensembl88.tsv?dl=1'
    path = str(my_data_pypath.join('protein_coding_genes_human_ensembl88.tsv'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_genome_annotation_file(my_data_pypath):
    _LOGGER.info('Downloading the genome annotation file...')
    url = 'https://www.dropbox.com/s/7oexg6k6aqukzi3/genome_annotations_human_ensembl88_test.gtf.gz?dl=1'
    path = str(my_data_pypath.join('genome_annotations_human_ensembl88_test.gtf.gz'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_barcode_counts_mapped_file(my_data_pypath):
    _LOGGER.info('Downloading the mapped barcode counts file...')
    url = 'https://www.dropbox.com/s/mq0v2dyskomjoon/SRR3879606_10k_barcode_counts_mapped.tsv?dl=1'
    path = str(my_data_pypath.join('barcode_counts_mapped.tsv'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_gene_expression_filtered_file(my_data_pypath):
    _LOGGER.info('Downloading the filtered gene expression file...')
    url = 'https://www.dropbox.com/s/1kixrh0cvgf5jqo/SRR3879606_gene_expression_filtered.tsv?dl=1'
    path = str(my_data_pypath.join('gene_expression_filtered.tsv'))
    misc.http_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_dummy_file(my_data_pypath):
    """An empty file."""
    path = str(my_data_pypath.join('dummy_file.txt'))
    pathlib.Path(path).touch()
    return path


@pytest.fixture(scope='session')
def my_dummy_dir(tmpdir_factory):
    """An empty directory."""
    dummy_dir = tmpdir_factory.mktemp('singlecell_dummy')
    return dummy_dir
