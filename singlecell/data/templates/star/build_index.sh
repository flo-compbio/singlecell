#!/bin/bash
# build STAR index

{% if echo is defined and echo %}
set -x  # echo on
{% endif %}

{% if use_docker is defined and use_docker %}
DOCKER_PREFIX="/host"
{% else %}
DOCKER_PREFIX=
{% endif %}

{% if use_docker is defined and use_docker %}
docker run -v "/:/host" quay.io/biocontainers/star:2.5.3a--0 \
{% endif %}
STAR \
    --runMode genomeGenerate \
    --runThreadN {{ num_threads }} \
    --genomeDir "${DOCKER_PREFIX}"'{{ output_dir }}' \
    --genomeFastaFiles "${DOCKER_PREFIX}"'{{ decompressed_genome_file }}' \
    {% if decompressed_genome_annotation_file is defined and 
            decompressed_genome_annotation_file is not none
        %}--sjdbGTFfile "${DOCKER_PREFIX}"'{{ decompressed_genome_annotation_file }}'
    {% endif %}
