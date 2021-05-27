# README

### To generate BAM files from synthetic variant genomes for testing variant callers

1. Create the virtual environment by running `./create.venv.sh` and activate it with `source mtb_amr-env/bin/activate`.

2. Simulate synthetic variant NGS reads by running `./01_simulate_variants.sh sample_name`

3. Preprocess and map reads by running `./02_preprocess_map_reads.sh sample_name`

4. Final output are `output/deduped/sample_name_deduped_mq10.bam` and `sample_name_deduped_mq10.bam.bai` for the next step of variant-calling. BAM files are also validated at this stage, and any errors can be found in `output/deduped/sample_name_deduped_mq10_validatebam.log`

*These two scripts should take ~10 mins to complete.*

### To call variants with all variant callers

1. Call variants by running `./03_call_variants.sh sample_name`

2. Variant callers currently available in the script: `gatk HaplotypeCaller`, `samtools/bcftools`, `freebayes`, `VarDict`, and `DiscoSnp`.

3. Final output is `output/called/sample_name_mq10_tool_normalized.vcf`. Output files from discosnp are contained in `output/called/discosnp/`.

*This script should take ~3 hrs to complete. (Calling variants with `VarDict` takes the most time.)*

### To check called variants against list of introduced variants

1. Check variants by running `./04_check_variants.sh sample_name`

2. Final output are `output/checked/sample_name_mq10_tool_normalized_fPOS.csv`, `output/checked/sample_name_mq10_tool_normalized_fNEG.csv`, and `output/checked/sample_name_mq10_tool_normalized_stats.csv`

*This script should take ~1 min to complete.*
