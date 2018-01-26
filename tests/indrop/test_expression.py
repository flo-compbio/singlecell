"""Tests for the `expression` module."""

import pytest

from genometools import misc

from singlecell.indrop import expression 

@pytest.mark.online
def test_quantify_expression(
        my_alignment_file, my_barcode_counts_mapped_file,
        my_chromosome_length_file, my_gene_file, my_genome_annotation_file,
        my_output_pypath):
    """Tests the gene expression quantification function."""

    #raw_output_file = \
    #        str(my_output_pypath.join('gene_expression_raw.mtx'))
    gene_expression_output_file = \
            str(my_output_pypath.join('gene_expression.mtx'))

    dense_gene_expression_output_file = \
            str(my_output_pypath.join('gene_expression.tsv'))

    transcript_expression_output_file = \
            str(my_output_pypath.join('transcript_expression.mtx'))

    num_cells = 10

    expression.quantify_expression(
        my_alignment_file, my_barcode_counts_mapped_file,
        my_chromosome_length_file, my_gene_file, my_genome_annotation_file,
        num_cells,
        gene_expression_output_file, transcript_expression_output_file,
        dense_gene_expression_output_file=dense_gene_expression_output_file,
        min_umi_qual=0,
        max_chroms=None)

    assert misc.get_file_md5sum(gene_expression_output_file) \
            == '7646339dea77247eedf330171bb87573'

    assert misc.get_file_md5sum(dense_gene_expression_output_file) \
            == '726320545f4755b3fbabecdffaa2dbed'

    assert misc.get_file_md5sum(transcript_expression_output_file) \
            == 'debf5fa9a227b939c25fb1929944d0f0'


@pytest.mark.skip(reason="no way of currently testing this")
def test_quantify_expression_UMI_quality(
        my_alignment_file, my_barcode_counts_mapped_file,
        my_chromosome_length_file, my_gene_file, my_genome_annotation_file,
        my_output_pypath):
    """Tests the gene expression quantification function."""

    #raw_output_file = \
    #        str(my_output_pypath.join('gene_expression_raw.mtx'))
    gene_expression_output_file = \
            str(my_output_pypath.join('gene_expression2.mtx'))

    dense_gene_expression_output_file = \
            str(my_output_pypath.join('gene_expression2.tsv'))

    transcript_expression_output_file = \
            str(my_output_pypath.join('transcript_expression2.mtx'))

    num_cells = 10

    expression.quantify_expression(
        my_alignment_file, my_barcode_counts_mapped_file,
        my_chromosome_length_file, my_gene_file, my_genome_annotation_file,
        num_cells,
        gene_expression_output_file, transcript_expression_output_file,
        dense_gene_expression_output_file=dense_gene_expression_output_file,
        min_umi_qual=10,
        max_chroms=None)

    assert misc.get_file_md5sum(gene_expression_output_file) \
            == '7646339dea77247eedf330171bb87573'

    assert misc.get_file_md5sum(dense_gene_expression_output_file) \
            == '726320545f4755b3fbabecdffaa2dbed'

    assert misc.get_file_md5sum(transcript_expression_output_file) \
            == 'debf5fa9a227b939c25fb1929944d0f0'
