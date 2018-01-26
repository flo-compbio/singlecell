#!/bin/bash
# align single-end reads with STAR

# => if use_docker=True, then we assume that all specified paths are absolute
# => if use_docker=True is not specified, then `STAR`` and `samtools``
#    executables must be in the PATH


{% if echo is defined and echo %}
set -x  # echo on
{% endif %}

READ_FILE='{{ read_file }}'
INDEX_DIR='{{ index_dir }}'
OUTPUT_PREFIX='{{ output_prefix }}'

# mkdir "$OUTPUT_PREFIX"

echo "Mapping reads for '${READ_FILE}' with STAR..."
echo "STAR index directory: '${INDEX_DIR}''"
echo "Output prefix: '${OUTPUT_PREFIX}'"

# FASTQ_FILE_NAME=${FASTQ_FILE##*/}
# FASTQ_FILE_NAME="${FASTQ_FILE##*/}"
# FASTQ_NAME="${FASTQ_FILE_NAME%%.*}"
{% if use_docker is defined and use_docker %}
DOCKER_PREFIX="/host"
{% else %}
DOCKER_PREFIX=
{% endif %}

{% if use_docker is defined and use_docker %}
docker run -v "/:/host" quay.io/biocontainers/star:2.5.3a--0 STAR \
{% else %}
STAR \
{% endif %} \
    --genomeDir "${DOCKER_PREFIX}{{ index_dir }}" \
    --readFilesIn "${DOCKER_PREFIX}{{ read_file }}" \
    --outFileNamePrefix "${DOCKER_PREFIX}{{ output_prefix }}" \
    --outMultimapperOrder "Random" \
    {% if num_threads is defined and num_threads is not none 
        %}--runThreadN {{ num_threads }} \
    {% endif
    %}{% if compressed is defined and compressed
        %}--readFilesCommand "zcat" \
    {% endif
    %}{% if out_mult_nmax is defined and out_mult_nmax is not none
        %}--outSAMmultNmax {{ out_mult_nmax }} \
    {% endif
    %}{% if out_sam_type is defined and out_sam_type is not none
        %}--outSAMtype {{ out_sam_type }} \
    {% endif
    %}{% if keep_unmapped is defined and keep_unmapped
        %}--outSAMunmapped Within \
    {% endif
    %}{% if seedSearchStartLmax is defined and seedSearchStartLmax is not none
        %}--seedSearchStartLmax {{ seedSearchStartLmax }} \
    {% endif
    %}{% if outFilterMismatchNmax is defined and outFilterMismatchNmax is not none
        %}--outFilterMismatchNmax {{ outFilterMismatchNmax }} \
    {% endif
    %}{% if outFilterMultimapNmax is defined and outFilterMultimapNmax is not none
        %}--outFilterMultimapNmax {{ outFilterMultimapNmax }} \
    {% endif
    %}{% if alignSJDBoverhangMin is defined and alignSJDBoverhangMin is not none
        %}--alignSJDBoverhangMin {{ alignSJDBoverhangMin }} \
    {% endif %}

{% if use_docker is defined and use_docker %}
docker run -v "/:/host" quay.io/biocontainers/samtools:1.4.1--0 samtools \
{% else %}
samtools \
{% endif %} \
    index "${DOCKER_PREFIX}${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"

{# --outSAMstrandField intronMotif \ #}