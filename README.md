# README

### Create virtual environment

1. Set up the virtual environment: `./create_venv.sh`.

2. Activate the virtual environment: `source mtb-pipeline-env/bin/activate`.

### Configure the SCons pipeline

1. Configure synthetic variant parameters by editing the values in `configs/variants_settings.conf`. These settings will be logged in `logs/variants/` once the pipeline is initiated.

2. Configure other parameters, including variant name, by editing the values in `configs/settings.conf`. 

### Run the SCons pipeline

1. Dry run: `scons -n`.

2. Run: `scons`

3. Debug run: `scons --debug=explain`

4. Variant callers included in the current pipeline: GATK HaplotypeCaller, bcftools, FreeBayes, DiscoSnp, DeepVariant, VarDict, Lancet.

5. A minimal reference genome is provided for testing in `data/NC_000962.3-small.fa`.

6. Note that DiscoSNP and Lancet require reference reads and BAM files to call variants. The following example reference files are provided for the minimal reference genome: 
	- `output/deduped/H37Rv-small_deduped_mq_H37Rv-small.bam`
	- `output/deduped/H37Rv-small_deduped_mq_H37Rv-small.bam.bai`
	- `output/reads/H37Rv-small.R1.trimmed.fq.gz`
	- `output/reads/H37Rv-small.R2.trimmed.fq.gz`

