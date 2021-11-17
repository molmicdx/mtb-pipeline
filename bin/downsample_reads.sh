#!/bin/bash

VARIANT=H37Rv_10-3DEL100X
DOWNSAMPLE=5/100
DS_VARIANT=H37Rv_10-3DELhead5X
READS_DIR=output/reads
LOG_DIR=logs/reads
BINDS=/mnt/disk2/molmicro,/mnt/disk15/molmicro,$PWD
SINGULARITY=/molmicro/common/singularity
SEQMAGICK_IMG=seqmagick-0.6.1.img

echo "[seqmagick]"
if [ ! -f "$LOG_DIR/${VARIANT}_seqmagick.log" ]; then
	singularity run -B $BINDS $SINGULARITY/$SEQMAGICK_IMG info $READS_DIR/${VARIANT}.R*.trimmed.fq.gz > $LOG_DIR/${VARIANT}_seqmagick.log
fi
echo "seqmagick done"

NUM_READS=$(cut -f 6 $LOG_DIR/${VARIANT}_seqmagick.log | tail -1)

let "DS_READS=$NUM_READS * $DOWNSAMPLE"

let "NUM_LINES=$DS_READS * 4"

echo "Downsampling $NUM_READS R1 reads to $DS_READS R1 reads..." 
gzip -cd $READS_DIR/${VARIANT}.R1.trimmed.fq.gz | head -n $NUM_LINES | gzip > $READS_DIR/${DS_VARIANT}.R1.trimmed.fq.gz
echo "Downsampling $NUM_READS R2 reads to $DS_READS R2 reads..." 
gzip -cd $READS_DIR/${VARIANT}.R2.trimmed.fq.gz | head -n $NUM_LINES | gzip > $READS_DIR/${DS_VARIANT}.R2.trimmed.fq.gz

echo "seq    num_reads" > $LOG_DIR/${DS_VARIANT}_downsampled.log
echo "$READS_DIR/${DS_VARIANT}.R1.trimmed.fq.gz    $NUM_READS" >> $LOG_DIR/${DS_VARIANT}_downsampled.log
echo "$READS_DIR/${DS_VARIANT}.R2.trimmed.fq.gz    $NUM_READS" >> $LOG_DIR/${DS_VARIANT}_downsampled.log
echo "Logged in $LOG_DIR/${DS_VARIANT}_downsampled.log"
