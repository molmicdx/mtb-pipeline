#!/bin/bash

# This script runs all the variant callers

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
SAMTOOLS=samtools-1.2.img
BCFTOOLS=bcftools-1.2.img
GATK=gatk-4.0.8.1.simg
FREEBAYES=quay.io/biocontainers/freebayes:1.3.5--py39hd2e4403_0
VARDICT=quay.io/biocontainers/vardict-java:1.8.2--0
VARDICT_SCRIPTS=bin/VarDict-1.8.2/bin
BEDTOOLS=bedtools-2.27.1-singularity-3.5.1.sif
DEDUPED_DIR=output/deduped
VC_DIR=output/called

# 10. Call variants with GATK HaplotypeCaller
echo "[GATK HaplotypeCaller] Calling variants..."
singularity exec -B $PWD $SINGULARITY/$GATK gatk HaplotypeCaller -ploidy 1 -R $REFERENCE_GENOME -I $DEDUPE_DIR/$1_deduped_matecig_mq10.bam -O $VC_DIR/$1_mq10_gatk.vcf > $VC_DIR/$1_mq10_gatk.log 2>&1
echo "Done"

# 11. Call variants with samtools/bcftools
echo "[samtools mpileup] Creating pileup file..."
singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools mpileup -m 3 -F 0.2 -u -f $REFERENCE_GENOME -d 100000 -A -B $DEDUPED_DIR/$1_deduped_matecig_mq10.bam -vo $VC_DIR/$1_mq10_samtools.vcf
echo "Done"

echo "[bcftools] Calling variants..."
singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools call -vmO v -o $VC_DIR/$1_mq10_bcftools.vcf $VC_DIR/$1_mq10_samtools.vcf
echo "Done"

# 12. Call variants with FreeBayes
echo "[FreeBayes] Calling variants..."
docker run -v `pwd`:`pwd` -w `pwd` -i -t --rm $FREEBAYES freebayes -f $REFERENCE_GENOME -p 1 --min-alternate-fraction 0.2 $DEDUPED_DIR/$1_deduped_matecig_mq10.bam > $VC_DIR/$1_mq10_freebayes.vcf
echo "Done"

# 13. Call variants with VarDict
echo "[bedtools] Creating BED file..."
singularity exec -B $PWD $SINGULARITY/$BEDTOOLS bedtools bamtobed -i $DEDUPED_DIR/$1_deduped_matecig_mq10.bam > $DEDUPED_DIR/$1_deduped_matecig_mq10.bed
echo "Done"

echo "[VarDict] Calling variants..."
docker run -v `pwd`:`pwd` -w `pwd` -i -t --rm $VARDICT vardict-java -G $REFERENCE_GENOME -f 0.2 -N $1 -b $DEDUPED_DIR/$1_deduped_matecig_mq10.bam -c 1 -S 2 -E 3 -g 4 $DEDUPED_DIR/$1_deduped_matecig_mq10.bed | $VARDICT_SCRIPTS/teststrandbias.R | $VARDICT_SCRIPTS/var2vcf_valid.pl -N $1 -E -f 0.2 > $VC_DIR/$1_mq10_vardict.vcf
echo "Done"

#14. Call variants with DiscoSnp
#echo "[DiscoSnp] Calling variants..."

#echo "Done"

#15. Normalize variant representation
echo "[GATK LeftAlignAndTrimVariants] Normalizing variant representations..."
singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_gatk.vcf -O $VC_DIR/$1_mq10_gatk_normalized.vcf > $VC_DIR/$1_mq10_gatk_normalized.log
singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_bcftools.vcf -O $VC_DIR/$1_mq10_bcftools_normalized.vcf > $VC_DIR/$1_mq10_bcftools_normalized.log
singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_freebayes.vcf -O $VC_DIR/$1_mq10_freebayes_normalized.vcf > $VC_DIR/$1_mq10_freebayes_normalized.log
singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_vardict.vcf -O $VC_DIR/$1_mq10_vardict_normalized.vcf > $VC_DIR/$1_mq10_vardict_normalized.log
echo "Done"

