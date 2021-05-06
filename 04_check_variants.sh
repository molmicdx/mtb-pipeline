#!/bin/bash

# This script checks variant calls against true mutations.

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
VARIANT_DIR=output/variants
VC_DIR=output/called
TO_CSV=bin/to_csv.py
CHECKER=bin/checker.py
CHECKED_DIR=output/checked

# Convert mutation list to csv

#python $TO_CSV $VARIANT_DIR/$1_normalized.vcf $1 --vcf

# Check true mutations

echo "Checking GATK variant calls..."
python $CHECKER $VC_DIR/$1_mq10_gatk_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv $CHECKED_DIR/$1_mq10_gatk_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_gatk_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_gatk_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking samtools variant calls..."
python $CHECKER $VC_DIR/$1_mq10_bcftools_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv $CHECKED_DIR/$1_mq10_bcftools_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_bcftools_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_bcftools_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking FreeBayes variant calls..."
python $CHECKER $VC_DIR/$1_mq10_freebayes_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv $CHECKED_DIR/$1_mq10_freebayes_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_freebayes_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_freebayes_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking VarDict variant calls..."
python $CHECKER $VC_DIR/$1_mq10_vardict_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv $CHECKED_DIR/$1_mq10_vardict_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_vardict_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_vardict_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking DiscoSnp variant calls..."
python $CHECKER $VC_DIR/discosnp/$1_discosnp-edit_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv $CHECKED_DIR/$1_mq10_discosnp-edit_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_discosnp-edit_normalized_fNEG.csv
mv $VC_DIR/discosnp/$1_discosnp-edit_normalized_stats.csv $CHECKED_DIR
echo "Done"

