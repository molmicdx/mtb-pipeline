#!/bin/bash
#This script creates reports on variant details

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_NAME=$1
DEDUPED_DIR=output/deduped
VC_DIR=output/called
IGV_REPORTS=output/igv_reports
INFO="AF DP MQ QD"

# Compress and index VCFs
echo "[bgzip | tabix] Compressing and indexing VCFs..."
for vcf in $VC_DIR/$1'_'*'_normalized.vcf'; do
	test -f $vcf && bgzip $PWD/$vcf;
	tabix -f $PWD/$vcf.gz;
done

if test -f $PWD/$VC_DIR/discosnp/$1_discosnp-edit_normalized_PASSsorted.vcf.gz; then
	echo "Already compressed"
else
	bgzip $PWD/$VC_DIR/discosnp/$1_discosnp-edit_normalized_PASSsorted.vcf
fi

if test -f $PWD/$VC_DIR/discosnp/$1_discosnp-edit_normalized_PASSsorted.vcf.gz.tbi; then
	echo "Already indexed"
else
	tabix -f $PWD/$VC_DIR/discosnp/$1_discosnp-edit_normalized_PASSsorted.vcf.gz
fi

if test -f $PWD/$VC_DIR/deepvariant/$1_mq10_deepvariant_normalized_PASS.vcf.gz; then
	echo "Already compressed"
else
	bgzip $PWD/$VC_DIR/deepvariant/$1_mq10_deepvariant_normalized_PASS.vcf
fi

if test -f $PWD/$VC_DIR/deepvariant/$1_mq10_deepvariant_normalized_PASS.vcf.gz.tbi; then 
	echo "Already indexed"
else
	tabix -f $PWD/$VC_DIR/deepvariant/$1_mq10_deepvariant_normalized_PASS.vcf.gz 
fi

echo "Done"

# Generate reports with igv-reports
echo "[igv-reports] Creating IGV report for GATK variant calls..."
create_report $VC_DIR/$1_mq10_gatk_normalized.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO --tracks $VC_DIR/$1_mq10_gatk_normalized.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_gatk_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

echo "[igv-reports] Creating IGV report for bcftools variant calls..."
create_report $VC_DIR/$1_mq10_bcftools_normalized.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO --tracks $VC_DIR/$1_mq10_bcftools_normalized.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_bcftools_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

echo "[igv-reports] Creating IGV report for FreeBayes variant calls..."
create_report $VC_DIR/$1_mq10_freebayes_normalized.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO --tracks $VC_DIR/$1_mq10_freebayes_normalized.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_freebayes_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

echo "[igv-reports] Creating IGV report for VarDict variant calls..."
create_report $VC_DIR/$1_mq10_vardict_normalized.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO --tracks $VC_DIR/$1_mq10_vardict_normalized.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_vardict_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

echo "[igv-reports] Creating IGV report for DiscoSnp variant calls..."
create_report $VC_DIR/discosnp/$1_discosnp-edit_normalized_PASSsorted.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO --tracks $VC_DIR/discosnp/$1_discosnp-edit_normalized_PASSsorted.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_discosnp-edit_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

echo "[igv-reports] Creating IGV report for DeepVariant variant calls..."
create_report $VC_DIR/deepvariant/$1_mq10_deepvariant_normalized_PASS.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO --tracks $VC_DIR/deepvariant/$1_mq10_deepvariant_normalized_PASS.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_deepvariant_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

echo "[igv-reports] Creating IGV report for delly variant calls..."
create_report $VC_DIR/$1_mq10_delly_normalized.vcf.gz $REFERENCE_GENOME --flanking 1000 --info-columns $INFO MAPQ SRMAPQ INSLEN --tracks $VC_DIR/$1_mq10_delly_normalized.vcf.gz $DEDUPED_DIR/$1_deduped_mq10.bam --output $IGV_REPORTS/$1_delly_report.html >> $IGV_REPORTS/igv_reports.log 2<&1
echo "Done"

