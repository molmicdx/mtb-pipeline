#!/bin/bash

# This script preprocesses and maps reads generated from sequencing of ATCC strain H37Ra by the Salipante lab, and maps it back to the H37Ra reference genome.
# The resulting indexed BAM files can then be used on variant callers.

REFERENCE_GENOME=data/CP000611.1.fasta
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
CUTADAPT=cutadapt-v1.15.img
BWA=bwa-0.7.13.img
SAMTOOLS=samtools-1.2.img
GATK=gatk-4.0.8.1.simg
INREADS_DIR=/home/local/AMC/ymseah/1089-710-517_atcc
READS_DIR=h37ra/output/reads
MAPPED_DIR=h37ra/output/mapped
DEDUPED_DIR=h37ra/output/deduped

# 5. Trim reads
echo "[cutadapt] Trimming reads..."
singularity exec -B $PWD $SINGULARITY/$CUTADAPT cutadapt --cores 8 -q 5 -b file:data/Illumina_adaptors_v2.fa -B file:data/Illumina_adaptors_v2.fa --minimum-length=20 -o $READS_DIR/$1.R1.trimmed.fq.gz -p $READS_DIR/$1.R2.trimmed.fq.gz $INREADS_DIR/$1_S4_L001_R1_001.fastq.gz $INREADS_DIR/$1_S4_L001_R2_001.fastq.gz > $READS_DIR/$1_cutadapt.log 2>&1
echo "Done"

#  6. Map reads
echo "[bwa mem] Mapping reads to reference genome..."
singularity exec -B $PWD $SINGULARITY/$BWA bwa mem $REFERENCE_GENOME $READS_DIR/$1.R1.trimmed.fq.gz $READS_DIR/$1.R2.trimmed.fq.gz -K 100000000 -R '@RG\tID:'$1'\tLB:LB_'$1'\tPL:illumina\tPU:TEST_RUN\tSM:'$1'' > $MAPPED_DIR/$1_trimmedRG_to_h37ra.sam
echo "Done"

# 7. Create sorted and indexed BAM
echo "[samtools] Creating sorted BAM files..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools view -bo $MAPPED_DIR/$1_trimmedRG_to_h37ra.bam $MAPPED_DIR/$1_trimmedRG_to_h37ra.sam
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools sort $MAPPED_DIR/$1_trimmedRG_to_h37ra.bam $MAPPED_DIR/$1_trimmedRG_to_h37ra-sorted
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools index $MAPPED_DIR/$1_trimmedRG_to_h37ra-sorted.bam
echo "Done"

# 8. Remove duplicate reads
echo "[gatk Picard MarkDuplicates] Removing duplicate reads..."
singularity exec -B $PWD $SINGULARITY/$GATK gatk MarkDuplicates -I $MAPPED_DIR/$1_trimmedRG_to_h37ra-sorted.bam -O $DEDUPED_DIR/$1_to_h37ra_deduped.bam -M $DEDUPED_DIR/$1_to_h37ra_deduped_metrics.txt --REMOVE_DUPLICATES TRUE > $DEDUPED_DIR/$1_to_h37ra_deduped.log 2>&1
echo "Done"

# 9. Filter based on mapping quality
echo "[samtools] Removing reads with low mapping quality..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools view $DEDUPED_DIR/$1_to_h37ra_deduped.bam -q 10 -bo $DEDUPED_DIR/$1_to_h37ra_deduped_mq10.bam
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools index $DEDUPED_DIR/$1_to_h37ra_deduped_mq10.bam
echo "Done"

# 10. Validate BAM
echo "[gatk ValidateSamFile] Validating BAM..."
singularity exec -B $PWD $SINGULARITY/$GATK/ gatk ValidateSamFile -I $DEDUPED_DIR/$1_to_h37ra_deduped_mq10.bam --MODE SUMMARY > $DEDUPED_DIR/$1_to_h37ra_deduped_mq10_validatebam.log 2<&1
echo "Done"


