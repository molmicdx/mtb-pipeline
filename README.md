# README

### To generate BAM files from synthetic variant genomes for testing variant callers

1. Create the virtual environment by running `./create.venv.sh`

2. Simulate synthetic variant NGS reads by running `./01_simulate_variants.sh sample_name`

3. Preprocess and map reads by running `./02_preprocess_map_reads.sh sample_name`

4. Final output is `output/deduped/sample_name_deduped_matecig_mq10.bam` and `sample_name_deduped_matecig_mq10.bam.bai`

### To call variants with all variant callers

1. Call variants by running `./03_call_variants.sh sample_name`

2. Final output is `output/called/sample_name_normalized.vcf`

