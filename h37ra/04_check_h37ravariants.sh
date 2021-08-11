#!/bin/bash

# This script checks variant calls on reads simulated from reference genomes. Expected number of variants is 0 since no mutations were synthetically introduced.

#REFERENCE_GENOME=data/CP000611.1.fasta
#REF=h37ra
REFERENCE_GENOME=data/CP016972.1.fasta
REF=CP016972
REF_NAME=$1
SINGULARITY=/molmicro/common/singularity
VARIANT_DIR=h37ra/output/variants
DEDUPED_DIR=h37ra/output/deduped
VC_DIR=h37ra/output/called
CHECKED_DIR=h37ra/output/checked
TO_CSV=bin/to_csv.py
TO_BED=bin/to_bed.py
FILTER_COV=bin/filter_low_cov.py
CHECKER=bin/checker.py
COV_LIMIT=10
BEDTOOLS=bedtools-2.27.1-singularity-3.5.1.sif

# Generate empty mutation files
#echo "CHROM,POS,REF,ALT,"$1 > $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv'
#touch $VARIANT_DIR/$1'_normalized.vcf.csv.bed'

#echo "[bedtools genomecov] Getting whole genome coverage..."
# https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
#singularity exec -B $PWD $SINGULARITY/$BEDTOOLS bedtools genomecov -ibam $DEDUPED_DIR/$1_to_$REF'_deduped_mq10.bam' -bga > $DEDUPED_DIR/$1_to_$REF'_deduped_mq10_genomecov.bed'
#echo "Done"

#echo "[bedtools intersect] Getting mutation coverage..."
# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
#singularity exec -B $PWD $SINGULARITY/$BEDTOOLS bedtools intersect -a $VARIANT_DIR/$1'_normalized.vcf.csv.bed' -b $DEDUPED_DIR/$1_to_$REF'_deduped_mq10_genomecov.bed' -wo | cut -f 1-3,7 > $VARIANT_DIR/$1_to_$REF'_normalized_genomecov_intersect.bed'
#echo "Done"

#echo "[filter_low_cov.py] Filter true mutations list..."
#python $FILTER_COV $VARIANT_DIR/$1_to_$REF'_normalized_genomecov_intersect.bed' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv' $COV_LIMIT
#echo "Done"

# Check true mutations
#echo "[checker.py] Checking GATK variant calls..."
#python $CHECKER $VC_DIR/$1_to_$REF'_mq10_gatk_normalized.vcf' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv_covfiltered.csv' $CHECKED_DIR/$1_to_$REF'_mq10_gatk_normalized_fPOS.csv' $CHECKED_DIR/$1_to_$REF'_mq10_gatk_normalized_fNEG.csv'
#echo "Done"

#echo "[checker.py] Checking samtools variant calls..."
#python $CHECKER $VC_DIR/$1_to_$REF'_mq10_bcftools_normalized.vcf' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv_covfiltered.csv' $CHECKED_DIR/$1_to_$REF'_mq10_bcftools_normalized_fPOS.csv' $CHECKED_DIR/$1_to_$REF'_mq10_bcftools_normalized_fNEG.csv'
#echo "Done"

#echo "[checker.py] Checking FreeBayes variant calls..."
#python $CHECKER $VC_DIR/$1_to_$REF'_mq10_freebayes_normalized.vcf' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv_covfiltered.csv' $CHECKED_DIR/$1_to_$REF'_mq10_freebayes_normalized_fPOS.csv' $CHECKED_DIR/$1_to_$REF'_mq10_freebayes_normalized_fNEG.csv'
#echo "Done"

echo "[checker.py] Checking VarDict variant calls..."
python $CHECKER $VC_DIR/$1_to_$REF'_mq10_vardict_normalized.vcf' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv_covfiltered.csv' $CHECKED_DIR/$1_to_$REF'_mq10_vardict_normalized_fPOS.csv' $CHECKED_DIR/$1_to_$REF'_mq10_vardict_normalized_fNEG.csv'
echo "Done"

echo "[checker.py] Checking DiscoSnp variant calls..."
python $CHECKER $VC_DIR/discosnp/$1_to_$REF'_discosnp-edit_normalized_PASSsorted.vcf' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv_covfiltered.csv' $CHECKED_DIR/$1_to_$REF'_discosnp-edit_normalized_PASSsorted_fPOS.csv' $CHECKED_DIR/$1_to_$REF'_discosnp-edit_normalized_PASSsorted_fNEG.csv'
echo "Done"

echo "[checker.py] Checking DeepVariant variant calls..."
python $CHECKER $VC_DIR/deepvariant/$1_to_$REF'_mq10_deepvariant_normalized_PASS.vcf' $VARIANT_DIR/$1_to_$REF'_normalized.vcf.csv_covfiltered.csv' $CHECKED_DIR/$1_to_$REF'_mq10_deepvariant_normalized_PASS_fPOS.csv' $CHECKED_DIR/$1_to_$REF'_mq10_deepvariant_normalized_PASS_fNEG.csv'
echo "Done"
