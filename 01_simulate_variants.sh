#!/bin/bash

# This script creates NGS reads from a simulated variant genome

VENV=mtb_amr-env
REFERENCE_GENOME=data/GCF_000195955.2_ASM19595v2_genomic.fna
VARIANT_SETTINGS=configs/variants_settings.conf
VARIANT_DIR=output/variants
VARIANT_NAME=$1
MNT_SINGULARITY=/molmicro,/molmicro/working/ymseah
SINGULARITY=/molmicro/common/singularity
GATK=gatk-4.0.8.1.simg
R1_QUAL_PROFILE=data/649-48_m_marinum_profileR1.txt
R2_QUAL_PROFILE=data/649-48_m_marinum_profileR2.txt
READ_LEN=150
AVG_READ_DEPTH=40
MEAN_FRAG_LEN=200
STD_DEV=10
READS_DIR=output/reads

# 1. Generate variant genome and list of introduced mutations
echo "[variants.py] Generating variant genome..."
cat $VARIANT_SETTINGS > $VARIANT_DIR/$1_variants_settings.conf
python bin/variants.py $REFERENCE_GENOME --settings $VARIANT_SETTINGS $VARIANT_DIR/$1.txt $VARIANT_DIR/$1.fa > $VARIANT_DIR/$1.log 2<&1
echo "Done"

# 2. Convert list of introduced mutations to VCF for normalization
echo "[to_vcf.py] Converting mutation list to VCF..."
python bin/to_vcf.py $VARIANT_DIR/$1.txt $1
echo "Done"

# 3. Normalize variant representation
echo "[gatk LeftAlignAndTrimVariants] Normalizing variant representation..."
singularity exec -B $MNT_SINGULARITY:/mnt $SINGULARITY/$GATK gatk LeftAlignAndTrimVariants -R /mnt/mtb_amr/$REFERENCE_GENOME -V /mnt/mtb_amr/$VARIANT_DIR/$1.txt.vcf -O /mnt/mtb_amr/$VARIANT_DIR/$1_normalized.vcf > $VARIANT_DIR/$1_normalized.log 2<&1
echo "Done"

# 4. Generate synthetic reads
echo "[art_illumina] Simulating reads..."
./$VENV/bin/art_illumina -1 $R1_QUAL_PROFILE -2 $R2_QUAL_PROFILE -p -sam -i $VARIANT_DIR/$1.fa -l $READ_LEN -f $AVG_READ_DEPTH -m $MEAN_FRAG_LEN -s $STD_DEV -o $READS_DIR/$1_R > $READS_DIR/$1_art_illumina.log 2<&1
echo "Compressing FASTQ files..."
gzip $READS_DIR/$1_R*.fq
echo "Done"

