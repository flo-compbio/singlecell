SingleCell
==========

| |pypi| |versions| |license|

..
    ===========  =================================================
    **latest**   |travis-latest| |codecov-latest| |docs-latest|
    **develop**  |travis-develop| |codecov-develop| |docs-develop|
    ===========  =================================================

SingleCell is a Python package for processing single-cell RNA-Seq data.

Requirements
------------

- Python 3 (tested with Python 3.5)
- STAR (tested with version 2.5.3a)
- samtools (tested with version 1.4.1)

The *STAR* and *samtools* executables must both be in the ``PATH``. To test
this, you can run the following commands, and check that they return the
respective version identifiers:

.. code-block:: bash
    
    $ STAR --version
    STAR_2.5.3a
    
    $ samtools --version
    samtools 1.4.1
    Using htslib 1.4.1
    Copyright (C) 2017 Genome Research Ltd.    

Installation
------------

.. code-block:: bash

    $ cd singlecell
    $ pip install -e .

Creating a STAR index (only once)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run the inDrop pipeline on your data, the first thing you need is a **STAR
genome index** for the species that your data is from. A STAR index consists of
a directory containing a bunch of files. For the human genome, the size of
these files totals about 25 GB. You only need to create an index once (per
species), which is then used by all future runs of the inDrop pipeline.

To generate an index, you need to download and decompress (using gunzip) the
genome (in FASTA file) and genome annotations (in GTF format) for the species
from the Ensembl FTP server. For example, for human:

.. code-block:: bash
    
    $ curl -O http://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    $ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    
    $ curl -O http://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz
    $ gunzip -c Homo_sapiens.GRCh38.88.gtf.gz > Homo_sapiens.GRCh38.88.gtf

For the genome annotation (GTF) file, you want to also keep the compressed
version, because this is the version used by the inDrop pipeline afterwards.

Now that you have those files ready, you can run the following:

.. code-block:: bash
    
    $ indrop_generate_star_index.py -g Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            -n Homo_sapiens.GRCh38.88.gtf \
            -od star_index_human -os build_star_index_human.sh \
            -ol build_star_index_human_log.txt \
            -t 16

This will output the STAR index in the directory "star_index_human"
(see ``-od`` parameter), and will use 16 threads in parallel (``-t``), making
the build process signficantly faster than if you were to run it
single-threaded.


Running the inDrop pipeline
---------------------------

To run the inDrop pipeline, you need to first create a configuration file (in
YAML format), which contains the locations (paths) of all the input files,
specifies an output directory, and sets a few parameters (e.g., how many cells
you want to include in the expression matrix). To generate a configuration file
template that you can then modify according to your setup, run the following:

.. code-block:: bash
    
    $ indrop_create_config_file.py -o my_configuration.yaml

After adjusting the parameters in the configuration file, you can check if
everything is configured correctly:

.. code-block:: bash
    
    $ indrop_check_pipeline.py -o my_configuration.yaml

If there are no errors, you can run the pipeline:

.. code-block:: bash
    
    $ indrop_pipeline.py -c my_configuration.yaml


.. |pypi| image:: https://img.shields.io/pypi/v/singlecell.svg
    :target: https://pypi.python.org/pypi/singlecell
    :alt: PyPI version

.. |versions| image:: https://img.shields.io/pypi/pyversions/singlecell.svg
    :target: https://pypi.python.org/pypi/singlecell
    :alt: Python versions supported

.. |license| image:: https://img.shields.io/pypi/l/singlecell.svg
    :target: https://pypi.python.org/pypi/singlecell
    :alt: License

.. |travis-latest| image:: https://travis-ci.org/flo-compbio/singlecell.svg?branch=master
    :alt: Build Status (master branch)
    :target: https://travis-ci.org/flo-compbio/singlecell

.. |travis-develop| image:: https://travis-ci.org/flo-compbio/singlecell.svg?branch=develop
    :alt: Build Status (develop branch)
    :target: https://travis-ci.org/flo-compbio/singlecell

.. |codecov-latest| image:: https://codecov.io/github/flo-compbio/singlecell/coverage.svg?branch=master
    :alt: Coverage (master branch)
    :target: https://codecov.io/github/flo-compbio/singlecell?branch=master

.. |codecov-develop| image:: https://codecov.io/github/flo-compbio/singlecell/coverage.svg?branch=develop
    :alt: Coverage (develop branch)
    :target: https://codecov.io/github/flo-compbio/singlecell?branch=develop

.. |docs-latest| image:: https://readthedocs.org/projects/singlecell/badge/?version=latest
    :alt: Documentation Status (master branch)
    :target: https://singlecell.readthedocs.org/en/latest

.. |docs-develop| image:: https://readthedocs.org/projects/singlecell/badge/?version=develop
    :alt: Documentation Status (develop branch)
    :target: https://singlecell.readthedocs.org/en/develop

.. _gtdocs: https://singlecell.readthedocs.org/en/latest/
