"""inDrop pipeline functions."""

import os
import logging
import time
import shutil
import logging

from pkg_resources import resource_filename
from plotly.offline import plot
import pandas as pd
import numpy as np

import snakemake

from jinja2 import Environment, PackageLoader, select_autoescape

_TEMPLATE_ENV = Environment(
    loader=PackageLoader('singlecell', os.path.join('data', 'indrop')),
    autoescape=select_autoescape(['html', 'xml']))


from genometools import misc
from genometools import ensembl
from genometools.expression import ExpMatrix

from .. import __version__
from .. import util
from .. import qc
from . import reads
from . import mapping
from . import aligned_reads
from . import barcodes
from . import expression
from . import config

_LOGGER = logging.getLogger(__name__)

def run_pipeline(config_file):
    
    conf, errors = config.read_config(config_file)

    pipeline = conf['pipeline']
    del conf['pipeline']

    output_dir = conf['output']['output_dir']

    template = _TEMPLATE_ENV.get_template(
        os.path.join('pipeline_template.sm'))

    misc.make_sure_dir_exists(output_dir)

    # write config file for pipeline
    config_file_name = 'config.yaml'
    pipeline_config_file = os.path.join(output_dir, config_file_name)
    config.write_config(conf, pipeline_config_file)

    # write snakemake file for pipeline
    snakemake_file = os.path.join(output_dir, 'Snakefile')
    with open(snakemake_file, 'w') as ofh:
        ofh.write(template.render(config_file=config_file_name))

    #with open(pipeline_config_file, 'w') as ofh:
    #    yaml.dump(conf, .write(conf, )

    num_threads = conf['parameters']['num_threads']
    snakemake.snakemake(snakemake_file, workdir=output_dir, cores=num_threads)
    

def run_pipeline_old(config_file):
    """inDrop pipeline."""

    t0 = time.time()

    conf, errors = config.read_config(config_file)

    input_ = conf['input']
    output = conf['output']
    params = conf['parameters']
    pipeline = conf['pipeline']

    output_dir = output['output_dir']

    barcode1_file = resource_filename(
        'singlecell', 'data/indrop/gel_barcode1_list.txt')
    barcode2_file = resource_filename(
        'singlecell', 'data/indrop/gel_barcode2_list.txt')

    if not util.is_empty_dir(output_dir):
        if output['allow_nonempty_output_dir']:
            _LOGGER.info('Note: Output directory exists and is not empty.')
        else:
            _LOGGER.error(
                'Output directory is not empty! Either specify an empty '
                '(or non-existent) output directory, or specify '
                '"allow_nonempty_output_dir: yes" in the configuration file.')
            return 1

    # create a timestamp for this run
    timestamp = time.strftime('%Y-%m-%d_%H-%M-%S')

    # create output directory, if necessary
    misc.make_sure_dir_exists(output_dir)

    # create results directory
    results_dir = os.path.join(output_dir, 'results')
    misc.make_sure_dir_exists(results_dir)

    # add file handler to the _LOGGER
    pipeline_log_file = os.path.join(results_dir, 'pipeline_log.txt')
    file_handler = logging.FileHandler(pipeline_log_file)
    log_fmt = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(log_fmt,log_datefmt)
    file_handler.setFormatter(formatter)
    _LOGGER.addHandler(file_handler)

    _LOGGER.info('This is the inDrop pipeline of SingleCell v%s', __version__)
    _LOGGER.info('Pipeline run timestamp: %s', timestamp)

    if params['use_docker']:
        _LOGGER.info('We\'re running using docker!')
    else:
        _LOGGER.info('We\'re running without using docker!')

    # create plot directory
    plot_dir = os.path.join(results_dir, 'qc_plots')
    misc.make_sure_dir_exists(plot_dir)

    # copy configuration file to results directory
    output_config_file = os.path.join(results_dir,
                                      'pipeline_config_%s.yaml' % timestamp)
    _LOGGER.info('Copying configuration file to "%s"', output_config_file)
    shutil.copyfile(config_file, output_config_file)
    
    ### process reads
    processed_read_dir = os.path.join(output_dir, 'processed_reads')
    misc.make_sure_dir_exists(processed_read_dir)
    process_read_file = os.path.join(
        processed_read_dir, 'processed_reads.fastq')
    process_count_file = os.path.join(results_dir, 'barcode_counts_reads.tsv')

    if pipeline['skip_read_processing']:
        _LOGGER.info('Skipping read processing step!')
    
    else:
        _LOGGER.info('Processing reads...')
        reads.process_reads(
            input_['barcode_read_file'], input_['mrna_read_file'],
            barcode1_file, barcode2_file,
            process_read_file, process_count_file,
            max_reads=params['max_reads'])
        _LOGGER.info('Finished processing reads.')

    
    ### mapping reads with STAR
    barcode_counts_mapped_file = os.path.join(results_dir,
                                              'barcode_counts_mapped.tsv')
    map_script_file = os.path.join(results_dir, 'map_with_star.sh')
    map_log_file = os.path.join(results_dir, 'mapping_log.txt')
    mapping_dir = os.path.join(output_dir, 'aligned_reads')

    alignment_file = os.path.join(mapping_dir, 'Aligned.sortedByCoord.out.bam')


    if pipeline['skip_mapping']:
        _LOGGER.info('Skipping read mapping step!')

    else:
        # mapping
        _LOGGER.info('Mapping reads with STAR...')
        star_params = conf['STAR']
        mapping.map_with_star(process_read_file, input_['star_index_dir'],
                              map_script_file, map_log_file,
                              mapping_dir,
                              num_threads=params['num_threads'],
                              compressed=False,
                              use_docker=params['use_docker'],
                              **star_params)
        _LOGGER.info('Finished mapping reads.')

        # count mapped reads for each barcode
        _LOGGER.info('Counting mapped reads for each barcode...')
        barcodes.count_mapped_reads(
            alignment_file,
            barcode1_file, barcode2_file,
            barcode_counts_mapped_file)
        _LOGGER.info('Finished counting mapped reads for each barcode.')

    ### generate intermediate files (chromosome lengths; protein-coding genes)
    chromlen_file = os.path.join(results_dir, 'chromosome_lengths.tsv')
    gene_file = os.path.join(results_dir, 'genes.tsv')
    if (not pipeline['skip_aligned_read_processing']) or \
            (not pipeline['skip_expression_quantification']):

        logger = logging.getLogger('genometools.ensembl')
        logger.setLevel(logging.ERROR)

        # generate file containing chromosome lengths
        _LOGGER.info('Extracting chromosome lengths...')
        chromlen = ensembl.get_chromosome_lengths(input_['genome_file'])
        chromlen.to_csv(chromlen_file, sep='\t', header=True)
        _LOGGER.info('Finished extracting chromosome lengths.')

        # generate file containing protein-coding genes
        _LOGGER.info('Extracting list of protein-coding genes from Ensembl GTF'
                     ' file...')
        protein_coding_genes = ensembl.get_protein_coding_genes(
            input_['genome_annotation_file'])
        _LOGGER.info('Finished extracting list of protein-coding genes.')

        if params['include_lincRNA_genes']:
            # extract lincRNA genes
            _LOGGER.info('Extracting list of lincRNA genes from Ensembl GTF'
                        ' file...')
            linc_rna_genes = ensembl.get_linc_rna_genes(
                input_['genome_annotation_file'])
            _LOGGER.info('Finished extracting list of lincRNA genes.')
            # exclude lincRNA whose gene name clashes with that of a
            # protein-coding gene
            sel = ~linc_rna_genes['name'].isin(
                set(protein_coding_genes['name']))
            linc_rna_genes = linc_rna_genes.loc[sel]
            genes = pd.concat([protein_coding_genes, linc_rna_genes])
        else:
            genes = protein_coding_genes
        
        genes.to_csv(gene_file, sep='\t', index=False)

    ### process aligned reads
    read_info_dir = os.path.join(output_dir, 'read_info')
    if not pipeline['skip_aligned_read_processing']:
        _LOGGER.info('Processing aligned reads...')
        misc.make_sure_dir_exists(read_info_dir)
        aligned_reads.process_aligned_reads(
            alignment_file, chromlen_file,
            gene_file, input_['genome_annotation_file'], read_info_dir,
            num_jobs=params['num_threads'])
        _LOGGER.info('Finished processing of aligned reads.')
    else:
        _LOGGER.info('Skipping processing of aligned reads!')

    ### quantify gene and transcript expression
    num_cells = params['num_cells']

    gene_expression_file = os.path.join(
        results_dir, 'gene_expression.mtx')
    transcript_expression_file = os.path.join(
        results_dir, 'transcript_expression.mtx')

    dense_gene_expression_file = None
    if output['generate_dense_expression_matrix']:
        dense_gene_expression_file = os.path.join(
            results_dir, 'gene_expression.tsv')

    if not pipeline['skip_expression_quantification']:
        _LOGGER.info('Quantifying expression for top %d cells...', num_cells)
        expression.quantify_expression(
            barcode_counts_mapped_file,
            chromlen_file,
            read_info_dir,
            gene_file, input_['genome_annotation_file'],
            num_cells,
            gene_expression_file, transcript_expression_file,
            min_umi_qual=params['min_umi_qual'],
            cell_prefix=output['cell_prefix'],
            dense_gene_expression_output_file=dense_gene_expression_file)
        _LOGGER.info('Finished expression quantification.')
    else:
        _LOGGER.info('Skipping expression quantification!')


    ### QC scripts
    if pipeline['skip_qc_plot_generation']:
        _LOGGER.info('Skipping the generation of QC plots!')

    else:
        _LOGGER.info('Generating QC plots...')

        _LOGGER.info('Plotting distribution of mapped reads per barcode...')
        mapped_count_histogram_file = \
                os.path.join(plot_dir, 'mapped_reads_histogram.html')
        barcodes.plot_read_count_distribution(
            barcode_counts_mapped_file, mapped_count_histogram_file)

        _LOGGER.info('Plotting distribution of transcripts per cell...')
        output_file = os.path.join(plot_dir, 'transcripts_per_cell.html')
        matrix = ExpMatrix.read_sparse(gene_expression_file)\
                .astype(np.float64)
        fig = qc.plot_cell_transcript_distribution(
            matrix, output['experiment_name'])
        plot(fig, filename=output_file, show_link=False, auto_open=False)

        _LOGGER.info('Plotting fraction of ribosomal and mitochondrial '
                     'gene expression...')
        output_file = os.path.join(plot_dir,
                                   'mito_ribo_expression.html')
        matrix = ExpMatrix.read_sparse(gene_expression_file)\
                .astype(np.float64)   # redundant
        fig = qc.plot_transcriptome_components(
            matrix, species=params['species'], name=output['experiment_name'],
            width=950, height=800, font_size=16, font_family='serif')
        plot(fig, filename=output_file, show_link=False, auto_open=False)

        _LOGGER.info('Plotting saturation...')
        output_file = os.path.join(plot_dir,
                                   'saturation.html')
        matrix = ExpMatrix.read_sparse(transcript_expression_file)\
                .astype(np.float64)
        fig = qc.plot_saturation(matrix)
        plot(fig, filename=output_file, show_link=False, auto_open=False)


    #_LOGGER.removeHandler(file_handler)
    t1 = time.time()
    t = t1 - t0
    _LOGGER.info('Pipeline run finished in %.1f s (%.1f min)!', t, t/60)
