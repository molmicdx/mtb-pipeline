#!/bin/bash

# This script checks variant calls against true mutations.

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
VARIANT_DIR=output/variants
DEDUPED_DIR=output/deduped
VC_DIR=output/called
CHECKED_DIR=output/checked
TO_CSV=bin/to_csv.py
TO_BED=bin/to_bed.py
FILTER_COV=bin/filter_low_cov.py
CHECKER=bin/checker.py
BEDTOOLS=bedtools-2.27.1-singularity-3.5.1.sif

echo "Creating CSV and BED files..."
# Convert mutation list to csv for checker.py
#python $TO_CSV $VARIANT_DIR/$1_normalized.vcf $1 --vcf

# Convert mutation list to bed for bedtools
python $TO_BED $VARIANT_DIR/$1_normalized.vcf.csv
echo "Done"

echo "[bedtools] Getting coverage depth at mutation positions..."
singularity exec -B $PWD $SINGULARITY/$BEDTOOLS bedtools genomecov -ibam $DEDUPED_DIR/$1_deduped_mq10.bam -bga > $DEDUPED_DIR/$1_deduped_mq10_genomecov.bed
singularity exec -B $PWD $SINGULARITY/$BEDTOOLS bedtools intersect -a $VARIANT_DIR/$1_normalized.vcf.csv.bed -b $DEDUPED_DIR/$1_deduped_mq10_genomecov.bed -wo | cut -f 1-3,7 > $VARIANT_DIR/$1_normalized_genomecov_intersect.bed
echo "Done"

echo "Filter true mutations list..."
python $FILTER_COV $VARIANT_DIR/$1_normalized_genomecov_intersect.bed $VARIANT_DIR/$1_normalized.vcf.csv 10
echo "Done"

# Check true mutations
echo "Checking GATK variant calls..."
python $CHECKER $VC_DIR/$1_mq10_gatk_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv_covfiltered.csv $CHECKED_DIR/$1_mq10_gatk_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_gatk_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_gatk_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking samtools variant calls..."
python $CHECKER $VC_DIR/$1_mq10_bcftools_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv_covfiltered.csv $CHECKED_DIR/$1_mq10_bcftools_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_bcftools_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_bcftools_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking FreeBayes variant calls..."
python $CHECKER $VC_DIR/$1_mq10_freebayes_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv_covfiltered.csv $CHECKED_DIR/$1_mq10_freebayes_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_freebayes_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_freebayes_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking VarDict variant calls..."
python $CHECKER $VC_DIR/$1_mq10_vardict_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv_covfiltered.csv $CHECKED_DIR/$1_mq10_vardict_normalized_fPOS.csv $CHECKED_DIR/$1_mq10_vardict_normalized_fNEG.csv
mv $VC_DIR/$1_mq10_vardict_normalized_stats.csv $CHECKED_DIR
echo "Done"

echo "Checking DiscoSnp variant calls..."
python $CHECKER $VC_DIR/discosnp/$1_discosnp-edit_normalized.vcf $VARIANT_DIR/$1_normalized.vcf.csv_covfiltered.csv $CHECKED_DIR/$1_discosnp-edit_normalized_fPOS.csv $CHECKED_DIR/$1_discosnp-edit_normalized_fNEG.csv
mv $VC_DIR/discosnp/$1_discosnp-edit_normalized_stats.csv $CHECKED_DIR
echo "Done"

