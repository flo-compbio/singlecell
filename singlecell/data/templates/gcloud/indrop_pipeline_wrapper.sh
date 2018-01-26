#!/bin/bash
# Wrapper for running inDrop pipeline on Google Cloud

# Author: Florian Wagner <florian.wagner@nyu.edu>

set -x  # echo on

### PARAMETERS

# input data parameters
SCRIPT_DIR='{{ script_dir }}'
BARCODE_READ_FILE='{{ barcode_read_file }}'
MRNA_READ_FILE='{{ mrna_read_file }}'
GENOME_FILE='{{ genome_file }}'
GENOME_ANNOTATION_FILE='{{ genome_annotation_file }}'
INDEX_DIR='{{ index_dir }}'  # star index
SINGLECELL_PKG_FILE='{{ singlecell_package_file }}'
SINGLECELL_PKG_VERSION='{{ singlecell_package_version }}'
PIPELINE_CONFIG_FILE='{{ pipeline_config_file }}'

# output parameters
OUTPUT_DIR='{{ output_dir }}'

# change to root's home directory
cd /root

# create subdirectories
mkdir scripts
mkdir download
mkdir data
mkdir output

### PREPARATION, PART 1: Install software

# download startup scripts
echo "Downloading startup scripts from Google Storage..."
gsutil -q cp -r "${SCRIPT_DIR}/*" scripts/
chmod -R a+x scripts/*

# install compiled crcmod module
echo "Installing compiled version of crcmod module..."
./scripts/install_crcmod.sh

# install docker
echo "Installing docker..."
./scripts/install_docker.sh

# install GCC
echo "Installing GCC..."
./scripts/install_gcc.sh

# install some missing packages
# (for HTSlib)
apt-get install -y zlib1g-dev libbz2-dev liblzma-dev

# install Python 3
echo "Installing Python 3..."
./scripts/install_python.sh

# install numpy
pip3 install numpy

# install cython
pip3 install cython

# install pandas
pip3 install pandas

# install genometools
pip3 install genometools


# download singlecell package and install

gsutil -q cp "$SINGLECELL_PKG_FILE" download/singlecell.tar.gz
gsutil -q cp "$PIPELINE_CONFIG_FILE" download/config.yaml

SINGLECELL_PKG_DIR="singlecell-${SINGLECELL_PKG_VERSION}"
tar -xzf download/singlecell.tar.gz -C download
pushd "download/${SINGLECELL_PKG_DIR}"
pip3 install -e .
popd


### PREPARATION, PART 2: Download all necessary data

# download input data
echo "Downloading the genome..."
gsutil -q cp "$GENOME_FILE" data/genome.fa.gz

echo "Downloading the genome annotations"
gsutil -q cp "$GENOME_ANNOTATION_FILE" data/genome_annotations.gtf.gz

{% if skip_read_processing %}
    echo "(we'll skip read processing)"
    {% if not skip_mapping %}
        echo "Downloading previously generated processed reads..."
        gsutil -q cp -r "${OUTPUT_DIR}/processed_reads" "output/"
    {% endif %}

{% else %}
    echo "Downloading the sequencing data..."
    gsutil -q cp "$BARCODE_READ_FILE" data/barcode_reads.fastq.gz
    gsutil -q cp "$MRNA_READ_FILE" data/mrna_reads.fastq.gz
{% endif %}

{% if skip_mapping %}
    echo "(we'll skip read mapping)"
    {% if not skip_aligned_read_processing %}
        echo "Downloading previously generated mapping results..."
        gsutil -q cp -r "${OUTPUT_DIR}/aligned_reads" "output/"
    {% endif %}

{% else %}
    echo "Downloading the STAR index..."
    gsutil -q cp -r "${INDEX_DIR}" data/
    INDEX_NAME=${INDEX_DIR##*/}
    mv "data/${INDEX_NAME}" "data/star_index"
{% endif %}

{% if skip_aligned_read_processing %}
    echo "(we'll skip aligned read processing)"
    {% if not skip_expression_quantification %}
        echo "Downloading previously generated read info data..."
        gsutil -q cp -r "${OUTPUT_DIR}/read_info" "output/"
    {% endif %}
{% endif %}

{% if skip_expression_quantification %}
    echo "(we'll skip expression quantification)"
    {% if not skip_qc_plot_generation %}
        echo "Downloading previously generated expression data..."
        gsutil -q cp -r "${OUTPUT_DIR}/results" "output/"
    {% endif %}
{% endif %}

# restart docker daemon
/etc/init.d/docker restart

# run the pipeline check
indrop_check_pipeline.py -c download/config.yaml

# run the pipeline
indrop_pipeline.py -c download/config.yaml

# copy the results
gsutil -q cp -r output/results "${OUTPUT_DIR}/"
gsutil -q cp -r output/read_info "${OUTPUT_DIR}/"
gsutil -q cp -r output/aligned_reads "${OUTPUT_DIR}/"
gsutil -q cp -r /var/log/syslog "${OUTPUT_DIR}/results/"

### RUN THE PIPELINE

# run the inDrop pipeline
#indrop_pipeline.py -c download/config.yaml

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
./scripts/delete_instance.sh
{% endif %}
