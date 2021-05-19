#!/bin/bash

# This script preprocesses and maps reads generated from sequencing of ATCC strain H37Ra by the Salipante lab, and maps it back to the H37Ra reference genome.
# The resulting indexed BAM files can then be used on variant callers.

REFERENCE_GENOME=data/CP000611.1.fasta
REF=h37ra
#REFERENCE_GENOME=data/CP016972.1.fasta
#REF=CP016972
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
CUTADAPT=cutadapt-v1.15.img
BWA=bwa-0.7.13.img
SAMTOOLS=samtools-1.2.img
GATK=gatk-4.0.8.1.simg
READS_DIR=h37ra/output/reads
MAPPED_DIR=h37ra/output/mapped
DEDUPED_DIR=h37ra/output/deduped

#  6. Map reads
echo "[bwa mem] Mapping reads to reference genome..."
singularity exec -B $PWD $SINGULARITY/$BWA bwa mem $REFERENCE_GENOME $READS_DIR/$1.R1.trimmed.fq.gz $READS_DIR/$1.R2.trimmed.fq.gz -K 100000000 -R '@RG\tID:'$1'\tLB:LB_'$1'\tPL:illumina\tPU:TEST_RUN\tSM:'$1 > $MAPPED_DIR/$1_trimmedRG_to_$REF.sam
echo "Done"

# 7. Create sorted and indexed BAM
echo "[samtools] Creating sorted BAM files..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools view -bo $MAPPED_DIR/$1_trimmedRG_to_$REF.bam $MAPPED_DIR/$1_trimmedRG_to_$REF.sam
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools sort $MAPPED_DIR/$1_trimmedRG_to_$REF.bam $MAPPED_DIR/$1_trimmedRG_to_$REF-sorted
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools index $MAPPED_DIR/$1_trimmedRG_to_$REF-sorted.bam
echo "Done"

# 8. Remove duplicate reads
echo "[gatk Picard MarkDuplicates] Removing duplicate reads..."
singularity exec -B $PWD $SINGULARITY/$GATK gatk MarkDuplicates -I $MAPPED_DIR/$1_trimmedRG_to_$REF-sorted.bam -O $DEDUPED_DIR/$1_to_$REF'_deduped.bam' -M $DEDUPED_DIR/$1_to_$REF'_deduped_metrics.txt' --REMOVE_DUPLICATES TRUE > $DEDUPED_DIR/$1_to_$REF'_deduped.log' 2>&1
echo "Done"

# 9. Filter based on mapping quality
echo "[samtools] Removing reads with low mapping quality..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools view $DEDUPED_DIR/$1_to_$REF'_deduped.bam' -q 10 -bo $DEDUPED_DIR/$1_to_$REF'_deduped_mq10.bam'
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools index $DEDUPED_DIR/$1_to_$REF'_deduped_mq10.bam'
echo "Done"

# 10. Validate BAM
echo "[gatk ValidateSamFile] Validating BAM..."
singularity exec -B $PWD $SINGULARITY/$GATK/ gatk ValidateSamFile -I $DEDUPED_DIR/$1_to_$REF'_deduped_mq10.bam' --MODE SUMMARY > $DEDUPED_DIR/$1_to_$REF'_deduped_mq10_validatebam.log' 2<&1
echo "Done"
