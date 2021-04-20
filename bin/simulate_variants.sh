#!/bin/bash

# This script creates NGS reads from a simulated variant genome

REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_SETTINGS=configs/variants_settings.conf
VARIANT_DIR=output/variants
VARIANT_NAME=$1
MNT_SINGULARITY=/molmicro,/molmicro/working/ymseah
SINGULARITY=/molmicro/common/singularity
GATK_IMG=gatk-4.0.8.1.simg
R1_QUAL_PROFILE=data/649-48_m_marinum_profileR1.txt
R2_QUAL_PROFILE=data/649-48_m_marinum_profileR2.txt

# 1. Generate variant genome and list of introduced mutations
echo "[variants.py] Generating variant genome..."
python bin/variants.py $REFERENCE_GENOME --settings $VARIANT_SETTINGS $VARIANT_DIR/$1.txt $VARIANT_DIR/$1.fa > $VARIANT_DIR/$1.log 2<&1
echo "Done"

# 2. Convert list of introduced mutations to VCF for normalization
echo "[to_vcf.py] Converting mutation list to VCF..."
python bin/to_vcf.py $VARIANT_DIR/$1.txt $1
echo "Done"

# 3. Normalize variant representation
echo "[gatk LeftAlignAndTrimVariants] Normalizing variant representation..."
singularity exec -B $MNT_SINGULARITY:/mnt $SINGULARITY/$GATK_IMG gatk LeftAlignAndTrimVariants -R /mnt/mtb_amr/$REFERENCE_GENOME -V /mnt/mtb_amr/$VARIANT_DIR/$1.txt.vcf -O /mnt/mtb_amr/$VARIANT_DIR/$1_normalized.vcf > $VARIANT_DIR/$1_normalized.log 2<&1
echo "Done"

# 4. Generate synthetic reads

