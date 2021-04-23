#!/bin/bash

# This script preprocesses and maps synthetic reads generated from 01_simulate_variants.sh.
# The resulting indexed BAM files can then be used on variant callers.

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
CUTADAPT=cutadapt-v1.15.img
BWA=bwa-0.7.13.img
SAMTOOLS=samtools-1.2.img
GATK=gatk-4.0.8.1.simg
READS_DIR=output/reads
MAPPED_DIR=output/mapped
DEDUPED_DIR=output/deduped

# 5. Trim reads
echo "[cutadapt] Trimming reads..."
singularity exec -B $PWD $SINGULARITY/$CUTADAPT cutadapt --cores 8 -q 5 -b file:data/Illumina_adaptors_v2.fa -B file:data/Illumina_adaptors_v2.fa --minimum-length=20 -o $READS_DIR/$1.R1.trimmed.fq.gz -p $READS_DIR/$1.R2.trimmed.fq.gz $READS_DIR/$1_R1.fq.gz $READS_DIR/$1_R2.fq.gz > $READS_DIR/$1_cutadapt.log 2>&1
echo "Done"

# 6. Map reads
echo "[bwa mem] Mapping reads to reference genome..."
singularity exec -B $PWD $SINGULARITY/$BWA bwa mem $REFERENCE_GENOME $READS_DIR/$1.R1.trimmed.fq.gz $READS_DIR/$1.R2.trimmed.fq.gz -K 100000000 -R '@RG\tID:'$1'\tLB:LB_'$1'\tPL:illumina\tPU:TEST_RUN\tSM:NC_000962.3' > $MAPPED_DIR/$1_trimmedRG.sam
echo "Done"

# 7. Create sorted BAM
echo "[samtools] Creating sorted BAM files..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools view -bo $MAPPED_DIR/$1_trimmedRG.bam $MAPPED_DIR/$1_trimmedRG.sam
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools sort $MAPPED_DIR/$1_trimmedRG.bam $MAPPED_DIR/$1_trimmedRG-sorted
echo "Done"

# 8. Remove duplicate reads
echo "[gatk Picard MarkDuplicatesWithMateCigar] Removing duplicate reads..."
singularity exec -B $PWD $SINGULARITY/$GATK gatk MarkDuplicatesWithMateCigar -I $MAPPED_DIR/$1_trimmedRG-sorted.bam -O $DEDUPED_DIR/$1_deduped_matecig.bam -M $DEDUPED_DIR/$1_deduped_matecig_metrics.txt --REMOVE_DUPLICATES TRUE > $DEDUPED_DIR/$1_deduped_matecig.log 2>&1
echo "Done"

# 9. Filter based on mapping quality
echo "[samtools] Removing reads with low mapping quality..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools view $DEDUPED_DIR/$1_deduped_matecig.bam -q 10 -bo $DEDUPED_DIR/$1_deduped_matecig_mq10.bam
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools index $DEDUPED_DIR/$1_deduped_matecig_mq10.bam
echo "Done"

