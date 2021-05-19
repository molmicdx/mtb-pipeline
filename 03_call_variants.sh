#!/bin/bash

# This script runs all the variant callers

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
REF_NAME=TEST-MTBREF
VARIANT_NAME=$1
SINGULARITY=/molmicro/common/singularity
SAMTOOLS=samtools-1.2.img
BCFTOOLS=bcftools-1.2.img
GATK=gatk-4.0.8.1.simg
FREEBAYES=quay.io/biocontainers/freebayes:1.3.5--py39hd2e4403_0
VARDICT=quay.io/biocontainers/vardict-java:1.8.2--0
VARDICT_SCRIPTS=bin/VarDict-1.8.2/bin
BEDTOOLS=bedtools-2.27.1-singularity-3.5.1.sif
DISCOSNP=discosnp_wrap-v2.2.10.simg
READS_DIR=output/reads
DEDUPED_DIR=output/deduped
VC_DIR=output/called

# 10. Call variants with GATK HaplotypeCaller
#echo "[GATK HaplotypeCaller] Calling variants..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk HaplotypeCaller -ploidy 1 -R $REFERENCE_GENOME -I $DEDUPED_DIR/$1_deduped_mq10.bam -O $VC_DIR/$1_mq10_gatk_no_gt.g.vcf -ERC GVCF > $VC_DIR/$1_mq10_gatk.log 2>&1
#echo "[GATK GenotypeGVCFs] Genotyping variants..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk GenotypeGVCFs -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_gatk_no_gt.g.vcf -O $VC_DIR/$1_mq10_gatk.vcf >> $VC_DIR/$1_mq10_gatk.log 2>&1
#echo "Done"

#echo "[GATK LeftAlignAndTrimVariants] Normalizing GATK variant representations..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_gatk.vcf -O $VC_DIR/$1_mq10_gatk_normalized.vcf > $VC_DIR/$1_mq10_gatk_normalized.log 2>&1
#echo "Done"

echo "[GATK VariantRecalibrator and ApplyVQSR] Filtering by VQSLOD..."
./bin/vqsr.sh $REFERENCE_GENOME $VC_DIR/$1_mq10_gatk_normalized.vcf
echo "Done"

# 11. Call variants with samtools/bcftools
#echo "[samtools mpileup] Creating pileup file..."
#singularity exec -B $PWD $SINGULARITY/$SAMTOOLS samtools mpileup -m 3 -F 0.2 -u -f $REFERENCE_GENOME -d 100000 -A -B $DEDUPED_DIR/$1_deduped_mq10.bam -vo $VC_DIR/$1_mq10_samtools.vcf
#echo "Done"

#echo "[bcftools] Calling variants..."
#singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools call -vmO v -o $VC_DIR/$1_mq10_bcftools.vcf $VC_DIR/$1_mq10_samtools.vcf
#echo "Done"

#echo "[GATK LeftAlignAndTrimVariants] Normalizing bcftools variant representations..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_bcftools.vcf -O $VC_DIR/$1_mq10_bcftools_normalized.vcf > $VC_DIR/$1_mq10_bcftools_normalized.log 2>&1
#echo "Done"

# 12. Call variants with FreeBayes
echo "[FreeBayes] Calling variants..."
#docker run -v `pwd`:`pwd` -w `pwd` -i -t --rm $FREEBAYES freebayes -f $REFERENCE_GENOME -p 1 --min-alternate-fraction 0.2 $DEDUPED_DIR/$1_deduped_mq10.bam > $VC_DIR/$1_mq10_freebayes.vcf
#docker run -v `pwd`:`pwd` -w `pwd` -i -t --rm $FREEBAYES freebayes -f $REFERENCE_GENOME -p 1 --min-alternate-fraction 0.2 --gvcf $DEDUPED_DIR/$1_deduped_mq10.bam > $VC_DIR/$1_mq10_freebayes.g.vcf
echo "Done"

echo "[GATK LeftAlignAndTrimVariants] Normalizing FreeBayes variant representations..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_freebayes.vcf -O $VC_DIR/$1_mq10_freebayes_normalized.vcf > $VC_DIR/$1_mq10_freebayes_normalized.log 2>&1
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_freebayes.g.vcf -O $VC_DIR/$1_mq10_freebayes_normalized.g.vcf >> $VC_DIR/$1_mq10_freebayes_normalized.log 2>&1
echo "Done"

# 13. Call variants with VarDict
#echo "[bedtools] Creating BED file..."
#singularity exec -B $PWD $SINGULARITY/$BEDTOOLS bedtools bamtobed -i $DEDUPED_DIR/$1_deduped_mq10.bam > $DEDUPED_DIR/$1_deduped_mq10.bed
#echo "Done"

#echo "[VarDict] Calling variants..."
#docker run -v `pwd`:`pwd` -w `pwd` -i -t --rm $VARDICT vardict-java -G $REFERENCE_GENOME -f 0.2 -N $1 -b $DEDUPED_DIR/$1_deduped_mq10.bam -c 1 -S 2 -E 3 -g 4 $DEDUPED_DIR/$1_deduped_mq10.bed | $VARDICT_SCRIPTS/teststrandbias.R | $VARDICT_SCRIPTS/var2vcf_valid.pl -N $1 -E -f 0.2 > $VC_DIR/$1_mq10_vardict.vcf
#echo "Done"

#echo "[GATK LeftAlignAndTrimVariants] Normalizing VarDict variant representations..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/$1_mq10_vardict.vcf -O $VC_DIR/$1_mq10_vardict_normalized.vcf > $VC_DIR/$1_mq10_vardict_normalized.log 2>&1
#echo "Done"

#14. Call variants with DiscoSnp
#echo "Creating file of files..."
#echo $PWD/$READS_DIR/$REF_NAME.R1.trimmed.fq.gz > $VC_DIR/discosnp/$REF_NAME'_fof.txt'
#echo $PWD/$READS_DIR/$REF_NAME.R2.trimmed.fq.gz >> $VC_DIR/discosnp/$REF_NAME'_fof.txt'

#echo $PWD/$READS_DIR/$1.R1.trimmed.fq.gz > $VC_DIR/discosnp/$1'_fof.txt'
#echo $PWD/$READS_DIR/$1.R2.trimmed.fq.gz >> $VC_DIR/discosnp/$1'_fof.txt'

#echo $REF_NAME'_fof.txt' > $VC_DIR/discosnp/$REF_NAME'_'$1'_fof.txt'
#echo $1'_fof.txt' >> $VC_DIR/discosnp/$REF_NAME'_'$1'_fof.txt'

#echo "[DiscoSnp] Calling variants..."
#singularity run --pwd $PWD -B $PWD $SINGULARITY/$DISCOSNP run_discoSnp++.sh -r $VC_DIR/discosnp/$REF_NAME'_'$1'_fof.txt' -P 6 -b 1 -k 31 -c auto -T -l -G $REFERENCE_GENOME -p $REF_NAME'_'$1 -u 12 > $VC_DIR/discosnp/$REF_NAME'_'$1'_discosnp.log' 2>&1

#mv $PWD/$REF_NAME'_'$1* $VC_DIR/discosnp/
#echo "Done"

#echo "[GATK LeftAlignAndTrimVariants] Normalizing DiscoSnp variant representations..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/discosnp/$REF_NAME'_'$1'_k_31_c_auto_D_100_P_6_b_1_coherent.vcf' -O $VC_DIR/discosnp/$REF_NAME'_'$1'_discosnp_normalized.vcf' > $VC_DIR/discosnp/$REF_NAME'_'$1'_discosnp_normalized.log' 2>&1
#echo "Done"

#echo "Formatting DiscoSnp VCF for checker.py/vcfpy..."
#sed -e 's/##SAMPLE/##sample/' -e 's/G1/'$REF_NAME'/' -e 's/G2/'$1'/' <$VC_DIR/discosnp/$REF_NAME'_'$1'_discosnp_normalized.vcf' >$VC_DIR/discosnp/$REF_NAME'_'$1'_discosnp-edit_normalized.vcf'

# Split VCFs into single-sample VCF (https://www.biostars.org/p/138694/#138783)
#for sample in $(zgrep -m 1 "^#CHROM" $VC_DIR'/discosnp/'$REF_NAME'_'$1'_discosnp-edit_normalized.vcf' | cut -f10-); do
#	singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools view -c1 -Ov -s $sample -o $VC_DIR'/discosnp/'$sample'_discosnp-edit_normalized.vcf'  $VC_DIR'/discosnp/'$REF_NAME'_'$1'_discosnp-edit_normalized.vcf'; done
#rm $VC_DIR'/discosnp/'$REF_NAME'_discosnp-edit_normalized.vcf'

# Filter mutations that PASS filter
#grep "#" $VC_DIR'/discosnp/'$1'_discosnp-edit_normalized.vcf' > $VC_DIR'/discosnp/'$1'_discosnp-edit_normalized_PASS.vcf'
#grep "$(printf '\t')PASS$(printf '\t')" $VC_DIR'/discosnp/'$1'_discosnp-edit_normalized.vcf' >> $VC_DIR'/discosnp/'$1'_discosnp-edit_normalized_PASS.vcf'

# Sort discosnp VCF
#singularity exec -B $PWD $SINGULARITY/$GATK gatk SortVcf -I $VC_DIR'/discosnp/'$1'_discosnp-edit_normalized_PASS.vcf' -O $VC_DIR'/discosnp/'$1'_discosnp-edit_normalized_PASSsorted.vcf'
#echo "Done"

# 15. Call variants with DeepVariant
#echo "[DeepVariant] Calling variants..."
#docker run -v $PWD:/input -v $PWD:/output --rm google/deepvariant:latest /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/input/$REFERENCE_GENOME --reads=/input/$DEDUPED_DIR/$1_deduped_mq10.bam --output_vcf=/output/$VC_DIR/deepvariant/$1_mq10_deepvariant.vcf --output_gvcf=/output/$VC_DIR/deepvariant/$1_mq10_deepvariant.gvcf --num_shards=12 --logging_dir=/output/$VC_DIR/deepvariant/
#echo "Done"

#echo "[GATK LeftAlignAndTrimVariants] Normalizing DeepVariant variant representations..."
#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/deepvariant/$1_mq10_deepvariant.vcf -O $VC_DIR/deepvariant/$1_mq10_deepvariant_normalized.vcf > $VC_DIR/deepvariant/$1_mq10_deepvariant_normalized.log 2>&1

#singularity exec -B $PWD $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R $REFERENCE_GENOME -V $VC_DIR/deepvariant/$1_mq10_deepvariant.gvcf -O $VC_DIR/deepvariant/$1_mq10_deepvariant_normalized.gvcf >> $VC_DIR/deepvariant/$1_mq10_deepvariant_normalized.log 2>&1
#echo "Done"

# Filter mutations that PASS filter
#echo "[grep] Get mutations that PASS..."
#grep "#" $VC_DIR'/deepvariant/'$1'_mq10_deepvariant_normalized.vcf' > $VC_DIR'/deepvariant/'$1'_mq10_deepvariant_normalized_PASS.vcf'
#grep "$(printf '\t')PASS$(printf '\t')" $VC_DIR'/deepvariant/'$1'_mq10_deepvariant_normalized.vcf' >> $VC_DIR'/deepvariant/'$1'_mq10_deepvariant_normalized_PASS.vcf'
#echo "Done"
