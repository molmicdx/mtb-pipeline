#!/bin/bash
# [Adapted from https://github.com/ksw9/mtb_variant_identification/scripts/vqsr.sh]
# Variant quality score recalibration (VQSR) with GATK. 
# Takes (1) reference genome and (2) VCF file.
# By default, this uses the variants with QUAL > mean QUAL as "truth set" for training. 

SINGULARITY=/molmicro/common/singularity
GATK=gatk-4.0.8.1.simg
BCFTOOLS=bcftools-1.2.img
VC_DIR=output/called/vqsr

# read from command line
ref=$1
vcf=$2

# Set up environment.
#source config.txt

# get basename
base=${vcf%.vcf*}
base=$(basename $base)

# get mean QUAL, excluding QUAL of invariant sites, to select variants for training
qual=$(singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools filter -i 'TYPE == "SNP"'  ${vcf}  | singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools query  -f '%QUAL\n' | grep -v inf | \
	awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')

# echo qual
echo 'mean qual of variant sites': $qual

# Create truth set by selecting only high qual variants (QUAL > mean) among variant sites.
singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools view --types snps ${vcf} | singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools filter -e 'QUAL == inf' | singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools filter -i " QUAL > $qual " > $VC_DIR/${base}_training_sites.vcf

# Index vcf.
bgzip -f $VC_DIR/${base}_training_sites.vcf > $VC_DIR/${base}_training_sites.vcf.gz
tabix -f -p vcf $VC_DIR/${base}_training_sites.vcf.gz

# Recalibrate variants. 

if singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output $VC_DIR/${base}.recal \
--tranches-file $VC_DIR/${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0 ; 
  then
 echo 'VQSR succeeded'

elif
# If failed, remove MQ which may not have sufficient variation.
singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-mode SNP \
--max-gaussians 2 \
--output $VC_DIR/${base}.recal \
--tranches-file $VC_DIR/${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0
  then
    echo 'VQSR w/o MQ'

elif
# If failed, remove MQRankSum which may not have sufficient variation.
singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output $VC_DIR/${base}.recal \
--tranches-file $VC_DIR/${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  
  then
    echo 'VQSR w/o MQRankSum'

elif
# If failed, remove ReadPosRankSum which may not have sufficient variation.
singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an MQRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output $VC_DIR/${base}.recal \
--tranches-file $VC_DIR/${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  
then 
  echo 'VQSR w/o ReadPosRankSum'

elif
      singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
	      -R ${ref} \
	      --variant ${vcf} \
	      --resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
	      -an DP \
	      -an QD \
	      -an FS \
	      -mode SNP \
	      --max-gaussians 2 \
	      --output $VC_DIR/${base}.recal \
	      --tranches-file $VC_DIR/${base}".tranches" \
	      -tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0
then
	echo 'VQSR w/o ReadPosRankSum or MQRankSum or MQ or SOR'
else 

# If failed, remove MQ, MQRankSum and ReadPosRankSum, which may not have sufficient variation.
singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
-an DP \
-an QD \
-an FS \
-an SOR \
-mode SNP \
--max-gaussians 2 \
--output $VC_DIR/${base}.recal \
--tranches-file $VC_DIR/${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0

echo 'VQSR w/o ReadPosRankSum or MQRankSum or MQ'

fi

# Set tranche filter level: this defines the sensitivity for the "truth variants."
ts_filter=99.0

# Apply VQSR.
if singularity exec -B $PWD $SINGULARITY/$GATK gatk ApplyVQSR \
-R ${ref} \
-mode SNP \
--variant ${vcf}  \
--recal-file $VC_DIR/${base}.recal  \
--tranches-file $VC_DIR/${base}".tranches"  \
--truth-sensitivity-filter-level ${ts_filter} \
--output $VC_DIR/${base}_vqsr.vcf
then
	echo "ApplyVQSR succeeded at ${ts_filter} ts_filter"
elif
	ts_filter=96.0
	singularity exec -B $PWD $SINGULARITY/$GATK gatk ApplyVQSR \
		-R ${ref} \
		-mode SNP \
		--variant ${vcf} \
		--recal-file $VC_DIR/${base}.recal \
		--tranches-file $VC_DIR/${base}".tranches" \
		--truth-sensitivity-filter-level ${ts_filter} \
		--output $VC_DIR/${base}_vqsr.vcf
then
	echo "ApplyVQSR succeeded at ${ts_filter} ts_filter"
elif
	ts_filter=93.0
	singularity exec -B $PWD $SINGULARITY/$GATK gatk ApplyVQSR \
		-R ${ref} \
		-mode SNP \
		--variant ${vcf} \
		--recal-file $VC_DIR/${base}.recal \
		--tranches-file $VC_DIR/${base}".tranches" \
		--truth-sensitivity-filter-level ${ts_filter} \
		--output $VC_DIR/${base}_vqsr.vcf
then
	echo "ApplyVQSR succeeded at ${ts_filter} ts_filter"
else
        ts_filter=90.0
	singularity exec -B $PWD $SINGULARITY/$GATK gatk ApplyVQSR \
		-R ${ref} \
		-mode SNP \
		--variant ${vcf} \
		--recal-file $VC_DIR/${base}.recal \
		--tranches-file $VC_DIR/${base}".tranches" \
		--truth-sensitivity-filter-level ${ts_filter} \
		--output $VC_DIR/${base}_vqsr.vcf

fi

# Tabix index and zip all VCF files for vcf-merge to work
bgzip -f -c $VC_DIR/${base}_vqsr.vcf >  $VC_DIR/${base}_vqsr.vcf.gz
tabix -f -p vcf $VC_DIR/${base}_vqsr.vcf.gz

# Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
echo 'VQSLOD contains inf/Nan at: ' $(singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools query -f '%INFO/VQSLOD\n' $VC_DIR/${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

# If VQSR does not run, set maxGuassians to 1. 
failedSites=$(singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools query  -f '%INFO/VQSLOD\n' $VC_DIR/${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)
if [ "$failedSites" -ne 0 ]; then
    echo 'rerunning with maxGaussians set to 1'
  
	singularity exec -B $PWD $SINGULARITY/$GATK gatk VariantRecalibrator \
	-R ${ref} \
	--variant ${vcf} \
	--resource truePos,known=false,training=true,truth=true,prior=15:$VC_DIR/${base}_training_sites.vcf.gz \
	-an DP \
	-an QD \
	-an MQRankSum \
	-an ReadPosRankSum \
	-an FS \
	-an SOR \
	-an MQ \
	-mode SNP \
	--max-gaussians 1 \
	--output $VC_DIR/${base}.recal \
	--tranches-file $VC_DIR/${base}".tranches"  \
	-tranche 100.0 -tranche 99.0 -tranche 99.9 -tranche 96.0 -tranche 93.0 -tranche 90.0 

	# define tranch filter level
	ts_filter=99.0

	# recalibration to SNPs 
	singularity exec -B $PWD $SINGULARITY/$GATK gatk ApplyVQSR \
	-R ${ref} \
	-mode SNP \
	--variant ${vcf}  \
	--recal-file $VC_DIR/${base}.recal  \
	--tranches-file $VC_DIR/${base}".tranches"  \
	--truth-sensitivity-filter-level ${ts_filter} \
	--output $VC_DIR/${base}_vqsr.vcf

	# tabix index and zip all VCF files for vcf-merge to work
	bgzip -f -c $VC_DIR/${base}_vqsr.vcf > $VC_DIR/${base}_vqsr.vcf.gz
	tabix -f -p vcf $VC_DIR/${base}_vqsr.vcf.gz
    
    # test if this works in calculating VQSLOD -- issue may be # gaussians is to high - use this to rewrite VQSR dir - change # gaussians if the VQSLOD is Nan/inf
    echo 'VQSLOD contains inf/Nan at: ' $(singularity exec -B $PWD $SINGULARITY/$BCFTOOLS bcftools query -f '%INFO/VQSLOD\n' $VC_DIR/${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

fi

# Remove intermediate files.
#rm ${base}.recal
#rm ${base}.recal.idx
#rm ${base}_training_sites.vcf.gz
#rm ${base}_training_sites.vcf.gz.tbi
#rm ${base}".tranches"
#rm ${base}_vqsr.vcf
