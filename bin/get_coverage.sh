#!/bin/bash
# Usage: bash get_coverage.sh <bam> <reference_accession>

BAM_FILE=$1
REF=$2
BINDS=/mnt/disk2/molmicro,/mnt/disk15/molmicro
SINGULARITY=/molmicro/common/singularity
BEDTOOLS=bedtools-2.27.1-singularity-3.5.1.sif

singularity exec -B $BINDS,$PWD $SINGULARITY/$BEDTOOLS bedtools genomecov -ibam $1 | grep $2 | awk '{ sum+=$2*$3/$4 } END { print sum }'
