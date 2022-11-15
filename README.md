# Identifying indels from WGS short reads of haploid genomes distinguishes variant-calling algorithms
## Abstract
Indel-containing reads are challenging to map to unique and correct genomic locations. We expand on previous efforts to validate bacterial variant-calling by developing a toolkit for introducing synthetic variants into the _Mycobacterium tuberculosis_ H37Rv reference genome, and measuring the precision and recall of seven variant callers in identifying SNPs, as well as small insertions and deletions. We used these _in silico_ altered genomes to test the performance of bcftools, DeepVariant, DiscoSNP, FreeBayes, GATK HaplotypeCaller, Lancet, and VarDict (Java implementation). These variant callers are widely used, have good documentation, and appear to be actively maintained within the last five years. Algorithms for identifying sequence variants from short reads can be broadly categorized into reference-based methods that rely on mapping reads onto a curated reference genome, and newer reference-free methods that construct deBruijn graphs of k-mers; some methods are a hybrid of both. All the variant callers we tested were reference-based methods, except for DiscoSNP.
## Pipeline
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

