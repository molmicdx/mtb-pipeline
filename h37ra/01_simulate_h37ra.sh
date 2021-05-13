#!/bin/bash

# This script simulates NGS reads from H37Ra (CP000611.1) reference genome

VENV=mtb_amr-env
REFERENCE_GENOME=data/CP000611.1.fasta
REF_NAME=TEST-CP000611
SINGULARITY=/molmicro/common/singularity
CUTADAPT=cutadapt-v1.15.img
READ_LEN=150
AVG_READ_DEPTH=40
MEAN_FRAG_LEN=200
STD_DEV=10
READS_DIR=h37ra/output/reads

# 4. Generate synthetic H37Ra reads from CP000611.1
echo "[art_illumina] Simulating reads..."
./$VENV/bin/art_illumina -p -sam -i $REFERENCE_GENOME -l $READ_LEN -f $AVG_READ_DEPTH -m $MEAN_FRAG_LEN -s $STD_DEV -o $READS_DIR/$REF_NAME'_R' > $READS_DIR/$REF_NAME'_art_illumina.log' 2<&1
echo "Compressing FASTQ files..."
gzip $READS_DIR/$REF_NAME'_R'*'.fq'
echo "Done"

# 5. Trim reads
echo "[cutadapt] Trimming reads..."
singularity exec -B $PWD $SINGULARITY/$CUTADAPT cutadapt --cores 8 -q 5 -b file:data/Illumina_adaptors_v2.fa -B file:data/Illumina_adaptors_v2.fa --minimum-length=20 -o $READS_DIR/$REF_NAME.R1.trimmed.fq.gz -p $READS_DIR/$REF_NAME.R2.trimmed.fq.gz $READS_DIR/$REF_NAME'_R1.fq.gz' $READS_DIR/$REF_NAME'_R2.fq.gz' > $READS_DIR/$REF_NAME'_cutadapt.log' 2>&1
echo "Done"
