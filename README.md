# README

### Create virtual environment

1. Set up the virtual environment: `./create_venv.sh`.

2. Activate the virtual environment: `source mtb_amr-env/bin/activate`.

### Configure the Scons pipeline

1. Configure synthetic variant parameters by editing the values in `configs/variants_settings.conf`. These settings will be logged in `logs/variants/` once the pipeline is initiated.

2. Configure other parameters, including variant name, by editing the values in `configs/settings.conf`. 

### Run the Scons pipeline

1. Dry run: `scons -n`.

2. Run: `scons`

3. Debug run: `scons --debug=explain`

4. Variant callers included in the current pipeline: GATK HaplotypeCaller, bcftools, FreeBayes, DiscoSnp, DeepVariant, VarDict, delly, Lancet.


### The pipeline in bash

This pipeline was initially put together with the following bash scripts:

#### To generate BAM files from synthetic variant genomes for testing variant callers

1. Simulate synthetic variant NGS reads by running `./01_simulate_variants.sh sample_name`

2. Preprocess and map reads by running `./02_preprocess_map_reads.sh sample_name`

3. Final output are `output/deduped/sample_name_deduped_mq10.bam` and `sample_name_deduped_mq10.bam.bai` for the next step of variant-calling. BAM files are also validated at this stage, and any errors can be found in `output/deduped/sample_name_deduped_mq10_validatebam.log`

*These two scripts should take ~10 mins to complete.*

#### To call variants with all variant callers

1. Call variants by running `./03_call_variants.sh sample_name`

2. Variant callers currently available in the script: `gatk HaplotypeCaller`, `samtools/bcftools`, `freebayes`, `VarDict`, `DiscoSnp`, and `DeepVariant`.

3. Final output is `output/called/sample_name_mq10_tool_normalized.vcf`. Output files from DiscoSnp and DeepVariant are contained in `output/called/discosnp/` and `output/called/deepvariant` respectively.

*This script should take ~3 hrs to complete. (Calling variants with `VarDict` takes the most time.)*

#### To check called variants against list of introduced variants

1. Check variants by running `./04_check_variants.sh sample_name`

2. Final output are `output/checked/sample_name_mq10_tool_normalized_fPOS.csv`, `output/checked/sample_name_mq10_tool_normalized_fNEG.csv`, and `output/checked/sample_name_mq10_tool_normalized_stats.csv`

*This script should take ~1 min to complete.*

#### To generate IGV reports

1. Generate IGV snapshot reports by running `./05_create_reports.sh sample_name`

2. Final HTML reports for each variant caller are in `output/igv_reports/`

*This script takes ~2 to ~45 mins to complete, depending on the number of variants*
