import os
import sys
from utils import parse_config

config = parse_config('configs/settings.conf')

# Ensure that a virtualenv is active before importing non-stdlib dependencies.
try:
    venv = os.environ.get('VIRTUAL_ENV')
except KeyError:
    sys.exit('an active virtualenv is required')

singularity = config.get('singularity', 'singularity')
gatk_img = config.get('singularity', 'gatk')
cutadapt_img = config.get('singularity', 'cutadapt')
bwa_img = config.get('singularity', 'bwa')
samtools_img = config.get('singularity', 'samtools')
bcftools_img = config.get('singularity', 'bcftools')
bedtools_img = config.get('singularity', 'bedtools')
discosnp_img = config.get('singularity', 'discosnp')
docker = config.get('docker', 'docker')
freebayes_img = config.get('docker', 'freebayes')
deepvariant_img = config.get('docker', 'deepvariant')
delly_img = config.get('singularity', 'delly')
vardict_img = config.get('docker', 'vardict')

from SCons.Script import (Environment, Variables, Help, Decider)

# check timestamps before calculating md5 checksums
Decider('MD5-timestamp')

# Define some PATH elements explicitly.
PATH=':'.join([
    'bin',
    os.path.join(venv, 'bin'),
    '/app/bin',  # provides R
    '/usr/local/bin', '/usr/bin', '/bin'])

vars = Variables()
vars.Add('out', '', 'output')
vars.Add('log', '', 'logs')
vars.Add('nproc', 'Number of concurrent processes', default=12)

# Provides access to options prior to instantiation of env object
# below; it's better to access variables through the env object.
# varargs = dict({opt.key: opt.default for opt in vars.options}, **vars.args)
# truevals = {True, 'yes', 'y', 'True', 'true', 't'}

# define boolean variables like
# varname = varargs['varname'] in truevals

# SHELLOPTS sets shell options to fail (including piped commands) with
# nonzero exit status; this requires bash.
env = Environment(
   # ENV = dict(os.environ, PATH=PATH, SHELLOPTS='errexit:pipefail'),
    ENV = dict(os.environ, PATH=PATH),
    variables = vars,
    SHELL = 'bash',
    cwd = os.getcwd(),
    venv_config = config.get('DEFAULT', 'virtualenv'),
    reference = config.get('DEFAULT', 'reference_genome'),
    ref_name = config.get('DEFAULT', 'reference_name'),
    accession = config.get('DEFAULT', 'reference_accession'),
    variant = config.get('variant_simulation', 'variant_name'),
    variants_config = config.get('variant_simulation', 'variants_config'),
    to_csv = config.get('variant_simulation', 'to_csv_script'),
    to_bed = config.get('variant_simulation', 'to_bed_script'),
    add_cov = config.get('variant_simulation', 'add_cov_script'),
    read1_q = config.get('variant_simulation', 'r1_qual_profile'),
    read2_q = config.get('variant_simulation', 'r2_qual_profile'),
    read_len = config.get('variant_simulation', 'read_length'),
    read_depth = config.get('variant_simulation', 'avg_read_depth'),
    frag_len = config.get('variant_simulation', 'avg_fragment_length'),
    sd = config.get('variant_simulation', 'std_dev'),
    min_read_q = config.get('reads_preprocessing', 'min_read_quality'),
    min_read_len = config.get('reads_preprocessing', 'min_read_length'),
    adaptors = config.get('reads_preprocessing', 'adaptors_fa'),
    cutadapt_cores = config.get('reads_preprocessing', 'cutadapt_cores'),
    variants_out = config.get('output', 'variants'),
    reads_out = config.get('output', 'reads'),
    mapped_out = config.get('output', 'mapped'),
    deduped_out = config.get('output', 'deduped'),
    called_out = config.get('output', 'called'),
    gvcf_out = config.get('output', 'gvcf'),
    gatk_out = config.get('variant_calling', 'gatk_output'),
    bcftools_out = config.get('variant_calling', 'bcftools_output'),
    ploidy_file = config.get('variant_calling', 'bcftools_ploidy_file'),
    bcftools_ann = config.get('variant_calling', 'bcftools_annotation'),
    freebayes_out = config.get('variant_calling', 'freebayes_output'),
    discosnp_out = config.get('variant_calling', 'discosnp_output'),
    deepvariant_out = config.get('variant_calling', 'deepvariant_output'),
    vardict_scripts = config.get('variant_calling', 'vardict_scripts'),
    vardict_out = config.get('variant_calling', 'vardict_output'),
    delly_out = config.get('variant_calling', 'delly_output'),
    lancet_out = config.get('variant_calling', 'lancet_output'),
    ploidy = config.get('variant_calling', 'ploidy'),
    allele_fraction = config.get('variant_calling', 'min_allele_fraction'),
    min_reads = config.get('variant_calling', 'min_reads'),
    max_reads = config.get('variant_calling', 'max_reads'),
    snp_per_bubble = config.get('discosnp_params', 'snp_per_bubble'),
    disco_mode = config.get('discosnp_params', 'disco_mode'),
    kmer_size = config.get('discosnp_params', 'kmer_size'),
    coverage = config.get('discosnp_params', 'coverage'),
    filter_low_complexity_bubbles = config.getboolean('discosnp_params', 'filter_low_complexity_bubbles'),
    max_threads = config.get('variant_calling', 'max_threads'),
    bwa_k = config.get('read_mapping', 'bwa_deterministic'),
    rg_pl = config.get('read_mapping', 'read_group_PL'),
    rg_pu = config.get('read_mapping', 'read_group_PU'),
    mapq = config.get('read_mapping', 'min_mapq'),
    bgz_out = config.get('output', 'bgzipped'),
    igv_out = config.get('output', 'igv_reports'),
    igv_info = config.get('report', 'igv_info'),
    igv_flank = config.get('report', 'igv_flanking_sites'),
    checked_out = config.get('output', 'checked'),
    check_call = config.get('report', 'checker_script'),
    min_read_depth = config.get('report', 'filter_read_depth'),
    mutation_to_flank = config.get('report', 'ignore_cov_on_mutation'), 
    num_flanking_bp = config.get('report', 'num_flanking_bp_to_use'),
    gatk = '{} exec -B $cwd {} gatk'.format(singularity, gatk_img),
    cutadapt = '{} exec -B $cwd {} cutadapt'.format(singularity, cutadapt_img),
    bwa = '{} exec -B $cwd {} bwa mem'.format(singularity, bwa_img),
    samtools = '{} exec -B $cwd {} samtools'.format(singularity, samtools_img),
    bcftools = '{} exec -B $cwd {}'.format(singularity, bcftools_img),
    bedtools = '{} exec -B $cwd {} bedtools'.format(singularity, bedtools_img),
    discosnp = '{} run --pwd $cwd -B $cwd {}'.format(singularity, discosnp_img),
    delly = '{} exec -B $cwd {} delly'.format(singularity, delly_img),
    freebayes = '{} run -v $cwd:$cwd -w $cwd -i -t --rm {} freebayes'.format(docker, freebayes_img),
    deepvariant = '{} run -v $cwd:/input -v $cwd:/output --rm {} /opt/deepvariant/bin/run_deepvariant'.format(docker, deepvariant_img),
    vardict = '{} run -v $cwd:$cwd -w $cwd -i -t --rm {} vardict-java'.format(docker, vardict_img)
)

# Help(vars.GenerateHelpText(env))

# ############### start inputs ################
# ############# Simulate Variants #############
simulated_variants_table, simulated_variants_fa = env.Command(
    target = ['$out/$variants_out/${variant}.txt', 
              '$out/$variants_out/${variant}.fa'],
    source = '$reference',
    action = ('python bin/variants.py --settings $variants_config $SOURCE $TARGETS > $log/$variants_out/${variant}.log 2>&1; ' 
              'cat $variants_config > $log/$variants_out/${variant}_variants_settings.conf')
)

# ############# Normalize Variant VCF ###############
simulated_variant_vcf = env.Command(
    target = '$out/$variants_out/${variant}.txt.vcf',
    source = simulated_variants_table,
    action = 'python bin/to_vcf.py $SOURCE $variant'
)

normalized_variant_vcf, normalized_log = env.Command(
    target = ['$out/$variants_out/${variant}_normalized.vcf', 
              '$log/$variants_out/${variant}_normalized.log'],
    source = simulated_variant_vcf,
    action = ('$gatk LeftAlignAndTrimVariants -R $reference -V $SOURCE -O ${TARGETS[0]} '
              '> ${TARGETS[-1]} 2>&1')
)

# ########### Simulate NGS Reads ##############
sim_R1, sim_R2, simreads_log = env.Command(
    target = ['$out/$reads_out/${variant}_R1.fq', 
              '$out/$reads_out/${variant}_R2.fq', 
              '$log/$reads_out/${variant}_art_illumina.log'],
    source = simulated_variants_fa,
    action = ('./${venv_config}/bin/art_illumina -1 $read1_q -2 $read2_q -p -sam '
              '-i $SOURCE -l $read_len -f $read_depth -m $frag_len -s $sd '
              '-o $out/$reads_out/${variant}_R > ${TARGETS[-1]} 2>&1')
)

gzsimR1, gzsimR2 = env.Command(
    target = ['$out/$reads_out/${variant}_R1.fq.gz',
              '$out/$reads_out/${variant}_R2.fq.gz'],
    source = [sim_R1, sim_R2],
    action = ('gzip $SOURCES --keep')
)

# ############## Trim NGS Reads ################
R1trimmed, R2trimmed, trimlog  = env.Command(
    target = ['$out/$reads_out/${variant}.R1.trimmed.fq.gz',
              '$out/$reads_out/${variant}.R2.trimmed.fq.gz',
              '$log/$reads_out/${variant}_cutadapt.log'],
    source = [gzsimR1, gzsimR2],
    action = ('$cutadapt --cores $cutadapt_cores -q $min_read_q -b file:$adaptors -B file:$adaptors '
              '--minimum-length $min_read_len -o ${TARGETS[0]} -p ${TARGETS[1]} $SOURCES '
              '> ${TARGETS[-1]} 2>&1')
)

# ################ Map Reads ###################

sam = env.Command(
    target = ['$out/$mapped_out/${variant}_trimmed.sam',
              '$log/$mapped_out/${variant}_trimmed_mapped.log'],
    source = ['$reference', R1trimmed, R2trimmed],
    action = ('$bwa $SOURCES -K $bwa_k '
              '-R \'@RG\\tID:${variant}\\tLB:LB_${variant}\\tPL:${rg_pl}\\tPU:${rg_pu}\\tSM:${variant}\' '
              '> ${TARGETS[0]} 2>${TARGETS[-1]}')
)

sorted_bam = env.Command(
    target = '$out/$mapped_out/${variant}_trimmed-sorted.bam',
    source = sam,
    action = '$samtools sort $SOURCE $out/$mapped_out/${variant}_trimmed-sorted'
)

# ############## Remove Duplicate Reads ###############
deduped, deduped_metrics, deduped_log = env.Command(
    target = ['$out/$deduped_out/${variant}_deduped.bam', 
              '$out/$deduped_out/${variant}_deduped_metrics.txt',
              '$log/$deduped_out/${variant}_deduped.log'],
    source = sorted_bam,
    action = ('$gatk MarkDuplicates -I $SOURCE -O ${TARGETS[0]} -M ${TARGETS[1]} --REMOVE_DUPLICATES TRUE '
              '> ${TARGETS[-1]} 2>&1')
)

# ################### Filter Reads, Get BAM Read Depths ####################

mq_filtered_bam, indexed_bam = env.Command(
    target = ['$out/$deduped_out/${variant}_deduped_mq.bam',
              '$out/$deduped_out/${variant}_deduped_mq.bam.bai'],
    source = deduped,
    action = ('$samtools view $SOURCE -q $mapq -bo ${TARGETS[0]}; '
              '$samtools index ${TARGETS[0]}')
)

genome_cov = env.Command(
    target = '$out/$variants_out/${variant}_deduped_mq10_genomecov.bed',
    source = mq_filtered_bam,
    action = '$bedtools genomecov -ibam $SOURCE -bga > $TARGET'
)

# ################## Validate BAM #####################

#bamlog = env.Command(
#    target = '$out/$deduped_out/${variant}_deduped_mq_validatebam.log',
#    source = mq_filtered_bam,
#    action = '$gatk ValidateSamFile -I $SOURCE --MODE SUMMARY > $TARGET 2>&1'
#)

# ############### end inputs ##################


# ################# Call Variants #####################

# ##################### GATK ##########################

gatk_gvcf, gatk_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${gatk_out}.g.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${gatk_out}_haplotypecaller.log'],
    source = ['$reference',
              mq_filtered_bam,
              indexed_bam],
    action = ('$gatk HaplotypeCaller -ploidy $ploidy -R ${SOURCES[0]} -I ${SOURCES[1]} '
              '-O ${TARGETS[0]} -ERC GVCF > ${TARGETS[-1]} 2>&1')
)

gatk_gt, gatk_gt_log = env.Command(
    target = ['$out/$called_out/$gatk_out/${variant}_${gatk_out}.vcf',
              '$log/$called_out/$gatk_out/${variant}_${gatk_out}_genotypegvcfs.log'],
    source = ['$reference',
              gatk_gvcf],
    action = ('$gatk GenotypeGVCFs -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

gatk_g_normalized, gatk_g_norm_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${gatk_out}_normalized.g.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${gatk_out}_normalized.g.log'],
    source = ['$reference',
              gatk_gvcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

gatk_normalized, gatk_norm_log = env.Command(
    target = ['$out/$called_out/$gatk_out/${variant}_${gatk_out}_normalized.vcf',
              '$log/$called_out/$gatk_out/${variant}_${gatk_out}_normalized.log'],
    source = ['$reference',
              gatk_gt],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

gatk_bgz = env.Command(
    target = '$out/$bgz_out/$gatk_out/${variant}_${gatk_out}_normalized.vcf.gz',
    source = gatk_normalized,
    action = 'bgzip < $SOURCE > $TARGET'
)

gatk_tbi, gatk_igv, gatk_igv_log = env.Command(
    target = ['$out/$bgz_out/$gatk_out/${variant}_${gatk_out}_normalized.vcf.gz.tbi',
              '$out/$igv_out/$gatk_out/${variant}_${gatk_out}_igv.html',
              '$log/$igv_out/$gatk_out/${variant}_${gatk_out}_igv.log'],
    source = [gatk_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

gatk_csv, gatk_bed = env.Command(
    target = ['$out/$called_out/$gatk_out/${variant}_${gatk_out}_normalized.vcf.csv',
              '$out/$called_out/$gatk_out/${variant}_${gatk_out}_normalized.vcf.csv.bed'],
    source = gatk_normalized,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)
gatk_cov_bed, gatk_cov_csv = env.Command(
    target = ['$out/$called_out/$gatk_out/${variant}_${gatk_out}_normalized_genomecov_intersect.bed',
              '$out/$called_out/$gatk_out/${variant}_${gatk_out}_normalized_cov.csv'],
    source = [gatk_csv,
              gatk_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ##################### bcftools #######################

bcftools_gvcf = env.Command(
    target = '$out/$called_out/$gvcf_out/${variant}_${bcftools_out}.g.vcf',
    source = ['$reference',
              mq_filtered_bam],
    action = ('$bcftools bcftools mpileup -Ou -f $SOURCES -d $max_reads -a $bcftools_ann 2>$log/$called_out/$gvcf_out/${variant}_${bcftools_out}.log | '
              'bcftools call -Ov -m --gvcf 1 --ploidy-file $ploidy_file -o $TARGET 2>>$log/$called_out/$gvcf_out/${variant}_${bcftools_out}.log')
)


bcftools_g_normalized, bcftools_g_norm_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${bcftools_out}_normalized.g.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${bcftools_out}_normalized.g.log'],
    source = ['$reference',
              bcftools_gvcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)


bcftools_vcf = env.Command(
    target = '$out/$called_out/$bcftools_out/${variant}_${bcftools_out}.vcf',
    source = ['$reference',
              mq_filtered_bam],
    action = ('$bcftools bcftools mpileup -Ou -f $SOURCES -d $max_reads -a $bcftools_ann 2>$log/$called_out/$bcftools_out/${variant}_${bcftools_out}.log | '
              'bcftools call -Ov -mv --ploidy-file $ploidy_file -o $TARGET 2>>$log/$called_out/$bcftools_out/${variant}_${bcftools_out}.log')
)

bcftools_normalized, bcftools_norm_log = env.Command(
    target = ['$out/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized.vcf',
              '$log/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized.log'],
    source = ['$reference',
              bcftools_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

bcftools_bgz = env.Command(
    target = '$out/$bgz_out/$bcftools_out/${variant}_${bcftools_out}_normalized.vcf.gz',
    source = bcftools_normalized,
    action = 'bgzip < $SOURCE > $TARGET'
)

bcftools_tbi, bcftools_igv, bcftools_igv_log = env.Command(
    target = ['$out/$bgz_out/$bcftools_out/${variant}_${bcftools_out}_normalized.vcf.gz.tbi',
              '$out/$igv_out/$bcftools_out/${variant}_${bcftools_out}_igv.html',
              '$log/$igv_out/$bcftools_out/${variant}_${bcftools_out}_igv.log'],
    source = [bcftools_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; '
              'create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

bcftools_csv, bcftools_bed = env.Command(
    target = ['$out/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized.vcf.csv',
              '$out/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized.vcf.csv.bed'],
    source = bcftools_normalized,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

bcftools_cov_bed, bcftools_cov_csv = env.Command(
    target = ['$out/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized_genomecov_intersect.bed',
              '$out/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized_cov.csv'],
    source = [bcftools_csv,
              bcftools_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ##################### FreeBayes  #######################

freebayes_gvcf = env.Command(
    target = '$out/$called_out/$gvcf_out/${variant}_${freebayes_out}.g.vcf',
    source = ['$reference',
              mq_filtered_bam],
    action = ('$freebayes -f ${SOURCES[0]} -p $ploidy --min-alternate-fraction $allele_fraction '
              '--gvcf ${SOURCES[1]} > $TARGET')
)

freebayes_vcf = env.Command(
    target = '$out/$called_out/$freebayes_out/${variant}_${freebayes_out}.vcf',
    source = ['$reference',
              mq_filtered_bam],
    action = ('$freebayes -f ${SOURCES[0]} -p $ploidy --min-alternate-fraction $allele_fraction '
              '${SOURCES[1]} > $TARGET')
)

freebayes_g_normalized, freebayes_g_norm_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${freebayes_out}_normalized.g.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${freebayes_out}_normalized.g.log'],
    source = ['$reference',
              freebayes_gvcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

freebayes_normalized, freebayes_norm_log = env.Command(
    target = ['$out/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized.vcf',
              '$log/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized.log'],
    source = ['$reference',
              freebayes_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

freebayes_bgz = env.Command(
    target = '$out/$bgz_out/$freebayes_out/${variant}_${freebayes_out}_normalized.vcf.gz',
    source = freebayes_normalized,
    action = 'bgzip < $SOURCE > $TARGET'
)

freebayes_tbi, freebayes_igv, freebayes_igv_log = env.Command(
    target = ['$out/$bgz_out/$freebayes_out/${variant}_${freebayes_out}_normalized.vcf.gz.tbi',
              '$out/$igv_out/$freebayes_out/${variant}_${freebayes_out}_igv.html',
              '$log/$igv_out/$freebayes_out/${variant}_${freebayes_out}_igv.log'],
    source = [freebayes_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; '
              'create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

freebayes_csv, freebayes_bed = env.Command(
    target = ['$out/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized.vcf.csv',
              '$out/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized.vcf.csv.bed'],
    source = freebayes_normalized,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

freebayes_cov_bed, freebayes_cov_csv = env.Command(
    target = ['$out/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized_genomecov_intersect.bed',
              '$out/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized_cov.csv'],
    source = [freebayes_csv,
              freebayes_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################# DeepVariant #################

deepvariant_gvcf, deepvariant_vcf, deepvariant_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${deepvariant_out}.g.vcf',
              '$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}.vcf',
              '$log/$called_out/$deepvariant_out/${variant}_${deepvariant_out}.log'],
    source = ['$reference',
              mq_filtered_bam],
    action = ('$deepvariant --model_type=WGS --ref=/input/${SOURCES[0]} --reads=/input/${SOURCES[1]} '
              '--output_gvcf=/output/${TARGETS[0]} --output_vcf=/output/${TARGETS[1]} --num_shards=$max_threads '
              '--logging_dir=/output/$log/$called_out/${deepvariant_out}/ > ${TARGETS[-1]} 2>&1')
)

deepvariant_g_normalized, deepvariant_g_norm_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${deepvariant_out}_normalized.g.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${deepvariant_out}_normalized.g.log'],
    source = ['$reference',
              deepvariant_gvcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

deepvariant_normalized, deepvariant_norm_log = env.Command(
    target = ['$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized.vcf',
              '$log/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized.log'],
    source = ['$reference',
              deepvariant_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

deepvariant_pass = env.Command(
    target = '$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_PASS.vcf',
    source = deepvariant_normalized,
    action = 'grep "#" $SOURCE > $TARGET; grep "$$(printf \'\\t\')PASS$$(printf \'\\t\')" $SOURCE >> $TARGET'
)

deepvariant_bgz = env.Command(
    target = '$out/$bgz_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_PASS.vcf.gz',
    source = deepvariant_pass,
    action = 'bgzip < $SOURCE > $TARGET'
)

deepvariant_tbi, deepvariant_igv, deepvariant_igv_log = env.Command(
    target = ['$out/$bgz_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_PASS.vcf.gz.tbi',
              '$out/$igv_out/$deepvariant_out/${variant}_${deepvariant_out}_igv.html',
              '$log/$igv_out/$deepvariant_out/${variant}_${deepvariant_out}_igv.log'],
    source = [deepvariant_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; '
              'create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

deepvariant_csv, deepvariant_bed = env.Command(
    target = ['$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized.vcf.csv',
              '$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized.vcf.csv.bed'],
    source = deepvariant_pass,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

deepvariant_cov_bed, deepvariant_cov_csv = env.Command(
    target = ['$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_genomecov_intersect.bed',
              '$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_cov.csv'],
    source = [deepvariant_csv,
              deepvariant_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################### Lancet ####################

lancet_vcf, lancet_log = env.Command(
    target = ['$out/$called_out/$lancet_out/${variant}_${ref_name}_${lancet_out}.vcf',
              '$log/$called_out/$lancet_out/${variant}_${ref_name}_${lancet_out}.log'],
    source = ['$reference',
              mq_filtered_bam,
              '$out/$deduped_out/${ref_name}_deduped_mq10.bam'],
    action = ('./bin/lancet --tumor ${SOURCES[1]} --normal ${SOURCES[2]} --ref ${SOURCES[0]} '
              '--reg $accession --min-vaf-tumor $allele_fraction --low-cov $min_read_depth '
              '--num-threads $max_threads --print-config-file > ${TARGETS[0]} 2>${TARGETS[-1]}')
)

lancet_config = env.Command(
    target = '$log/$called_out/$lancet_out/${variant}_${ref_name}_${lancet_out}_config.txt',
    source = None,
    action = 'mv ${cwd}/config.txt $TARGET'
)

lancet_normalized, lancet_norm_log = env.Command(
    target = ['$out/$called_out/$lancet_out/${variant}_${ref_name}_${lancet_out}_normalized.vcf',
              '$log/$called_out/$lancet_out/${variant}_${ref_name}_${lancet_out}_normalized.log'],
    source = ['$reference',
              lancet_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

lancet_final_vcf, lancet_sort_log = env.Command(
    target = ['$out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted.vcf',
              '$log/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted.log'],
    source = lancet_normalized,
    action = ('for sample in $$(zgrep -m 1 "^#CHROM" $SOURCE | cut -f10-); do '
              '    $bcftools bcftools view -c 1 -Ov -s $$sample -o $out/$called_out/$lancet_out/$$sample\'_${lancet_out}_normalized.vcf\' $SOURCE; done; '
              '$gatk SortVcf -I $out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized.vcf -O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1; '
              'rm $out/$called_out/$lancet_out/${ref_name}_${lancet_out}_normalized.vcf')
)


lancet_bgz = env.Command(
    target = '$out/$bgz_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted.vcf.gz',
    source = lancet_final_vcf,
    action = 'bgzip < $SOURCE > $TARGET'
)

lancet_tbi, lancet_igv, lancet_igv_log = env.Command(
    target = ['$out/$bgz_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted.vcf.gz.tbi',
              '$out/$igv_out/$lancet_out/${variant}_${lancet_out}_igv.html',
              '$log/$igv_out/$lancet_out/${variant}_${lancet_out}_igv.log'],
    source = [lancet_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; '
              'create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

lancet_csv, lancet_bed = env.Command(
    target = ['$out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted.vcf.csv',
              '$out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted.vcf.csv.bed'],
    source = lancet_final_vcf,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

lancet_cov_bed, lancet_cov_csv = env.Command(
    target = ['$out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted_genomecov_intersect.bed',
              '$out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_sorted_cov.csv'],
    source = [lancet_csv,
              lancet_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################## DiscoSnp ###################

ref_fof = env.Command(
    target = '$out/$called_out/$discosnp_out/${ref_name}_fof.txt',
    source = [R1trimmed, R2trimmed],
    action = ('echo ${cwd}/${out}/${reads_out}/${ref_name}.R1.trimmed.fq.gz > $TARGET; '
              'echo ${cwd}/${out}/${reads_out}/${ref_name}.R2.trimmed.fq.gz >> $TARGET')
)

variant_fof = env.Command(
    target = '$out/$called_out/$discosnp_out/${variant}_fof.txt',
    source = None,
    action = ('echo ${cwd}/${out}/${reads_out}/${variant}.R1.trimmed.fq.gz > $TARGET; '
              'echo ${cwd}/${out}/${reads_out}/${variant}.R2.trimmed.fq.gz >> $TARGET')
)

fof = env.Command(
    target = '$out/$called_out/$discosnp_out/${ref_name}_${variant}_fof.txt',
    source = [ref_fof,
              variant_fof],
    action = ('echo ${ref_name}_fof.txt > $TARGET; '
              'echo ${variant}_fof.txt >> $TARGET')
)

discosnp_vcf, discosnp_log = env.Command(
    target = ['$out/$called_out/$discosnp_out/${ref_name}_${variant}_k_${kmer_size}_c_${coverage}_D_100_P_${snp_per_bubble}_b_${disco_mode}_coherent.vcf',
              '$log/$called_out/${ref_name}_${variant}_discosnp.log'],    
    source = ['$reference',
              fof],
    action = ('$discosnp $out -r ../${SOURCES[1]} -P $snp_per_bubble '
              '-b $disco_mode -k $kmer_size -c $coverage -T -l '
              '-G ../${SOURCES[0]} -p ${ref_name}_${variant} -u $max_threads > ${TARGETS[-1]} 2>&1; '
              'mv $out/${ref_name}_${variant}* $out/$called_out/$discosnp_out/; '
	      'sed -i \'s/INDEL_.*_path_[0-9]*/${accession}/g\' ${TARGETS[0]}') #temp fix for inexplicable VCF entry; may falsely inflate false positive calls
)

discosnp_normalized, discosnp_norm_log = env.Command(
    target = ['$out/$called_out/$discosnp_out/${ref_name}_${variant}_discosnp_normalized.vcf',
              '$log/$called_out/$discosnp_out/${ref_name}_${variant}_discosnp_normalized.log'],
    source = ['$reference',
              discosnp_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

discosnp_formatted = env.Command(
    target = '$out/$called_out/$discosnp_out/${ref_name}_${variant}_discosnp-edit_normalized.vcf',
    source = discosnp_normalized,
    action = ('sed -e \'s/##SAMPLE/##sample/\' -e \'s/G1/${ref_name}/\' -e \'s/G2/${variant}/\' '
              '< $SOURCE > $TARGET')
)

discosnp_final_vcf = env.Command(
    target = '$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized.vcf',
    source = discosnp_formatted,
    action = ('for sample in $$(zgrep -m 1 "^#CHROM" $SOURCE | cut -f10-); do '
              '    $bcftools bcftools view -c 1 -Ov -s $$sample -o $out/$called_out/$discosnp_out/$$sample\'_discosnp-edit_normalized.vcf\' $SOURCE; done; '
              'rm $out/$called_out/$discosnp_out/${ref_name}_discosnp-edit_normalized.vcf')
)

discosnp_pass = env.Command(
    target = '$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASS.vcf',
    source = discosnp_final_vcf,
    action = 'grep "#" $SOURCE > $TARGET; grep "$$(printf \'\\t\')PASS$$(printf \'\\t\')" $SOURCE >> $TARGET'
)

discosnp_pass_sorted, discosnp_ps_log = env.Command(
    target = ['$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted.vcf',
              '$log/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted.log'],
    source = discosnp_pass,
    action = '$gatk SortVcf -I $SOURCE -O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1 '
)

discosnp_bgz = env.Command(
    target = '$out/$bgz_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted.vcf.gz',
    source = discosnp_pass_sorted,
    action = 'bgzip < $SOURCE > $TARGET'
)

discosnp_tbi, discosnp_igv, discosnp_igv_log = env.Command(
    target = ['$out/$bgz_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted.vcf.gz.tbi',
              '$out/$igv_out/$discosnp_out/${variant}_${discosnp_out}_igv.html',
              '$log/$igv_out/$discosnp_out/${variant}_${discosnp_out}_igv.log'],
    source = [discosnp_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; '
              'create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

discosnp_csv, discosnp_bed = env.Command(
    target = ['$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted.vcf.csv',
              '$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted.vcf.csv.bed'],
    source = discosnp_pass_sorted,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

discosnp_cov_bed, discosnp_cov_csv = env.Command(
    target = ['$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted_genomecov_intersect.bed',
              '$out/$called_out/$discosnp_out/${variant}_discosnp-edit_normalized_PASSsorted_cov.csv'],
    source = [discosnp_csv,
              discosnp_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################### VarDict ####################

bedfile = env.Command(
    target = '$out/$deduped_out/${variant}_deduped_mq10.bed',
    source = mq_filtered_bam,
    action = '$bedtools bamtobed -i $SOURCE > $TARGET'
)

vardict_vcf = env.Command(
    target = '$out/$called_out/$vardict_out/${variant}_vardict.vcf',
    source = ['$reference',
              mq_filtered_bam,
              bedfile],
    action = ('$vardict -G ${SOURCES[0]} -f $allele_fraction -N $variant -b ${SOURCES[1]} '
              '-c 1 -S 2 -E 3 -g 4 ${SOURCES[2]} | ${vardict_scripts}/teststrandbias.R '
              '| ${vardict_scripts}/var2vcf_valid.pl -N $variant -E -f $allele_fraction '
              '> $TARGET')
)

vardict_normalized, vardict_norm_log = env.Command(
    target = ['$out/$called_out/$vardict_out/${variant}_${vardict_out}_normalized.vcf',
              '$log/$called_out/$vardict_out/${variant}_${vardict_out}_normalized.log'],
    source = ['$reference',
              vardict_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

vardict_bgz = env.Command(
    target = '$out/$bgz_out/$vardict_out/${variant}_${vardict_out}_normalized.vcf.gz',
    source = vardict_normalized,
    action = 'bgzip < $SOURCE > $TARGET'
)

vardict_tbi, vardict_igv, vardict_igv_log = env.Command(
    target = ['$out/$bgz_out/$vardict_out/${variant}_${vardict_out}_normalized.vcf.gz.tbi',
              '$out/$igv_out/$vardict_out/${variant}_${vardict_out}_igv.html',
              '$log/$igv_out/$vardict_out/${variant}_${vardict_out}_igv.log'],
    source = [vardict_bgz,
              '$reference',
              mq_filtered_bam],
    action = ('tabix -f ${SOURCES[0]}; '
              'create_report ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]} > ${TARGETS[-1]} 2>&1')
)

vardict_csv, vardict_bed = env.Command(
    target = ['$out/$called_out/$vardict_out/${variant}_${vardict_out}_normalized.vcf.csv',
              '$out/$called_out/$vardict_out/${variant}_${vardict_out}_normalized.vcf.csv.bed'],
    source = vardict_normalized,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

vardict_cov_bed, vardict_cov_csv = env.Command(
    target = ['$out/$called_out/$vardict_out/${variant}_${vardict_out}_normalized_genomecov_intersect.bed',
              '$out/$called_out/$vardict_out/${variant}_${vardict_out}_normalized_cov.csv'],
    source = [vardict_csv,
              vardict_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ############### Get Variant Coverage ################

variant_csv, variant_bed = env.Command(
    target = ['$out/$variants_out/${variant}_normalized.vcf.csv',
              '$out/$variants_out/${variant}_normalized.vcf.csv.bed'],
    source = normalized_variant_vcf,
    action = ('python $to_csv $SOURCE ${TARGETS[0]}; '
              'python $to_bed $TARGETS --split_mut $mutation_to_flank --bp $num_flanking_bp')
)

variant_cov_bed, variant_cov_csv = env.Command(
    target = ['$out/$variants_out/${variant}_normalized_genomecov_intersect.bed',
              '$out/$variants_out/${variant}_normalized_cov.csv'],
    source = [variant_csv,
              variant_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo > ${TARGETS[0]}; '
              'python $add_cov ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################### Check Calls ######################

gatk_cov_filtered = env.Command(
    target = '$out/$called_out/$gatk_out/${variant}_${gatk_out}_normalized_dp${min_read_depth}.vcf',
    source = gatk_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

gatk_calls = env.Command(
    target = '$out/$checked_out/$gatk_out/${variant}_${gatk_out}_normalized_dp${min_read_depth}_checked.csv',
    source = [gatk_cov_filtered,
              variant_cov_csv,
              gatk_cov_csv],
    action = 'python $check_call $gatk_out $variant $SOURCES $TARGET'
)


bcftools_cov_filtered = env.Command(
    target = '$out/$called_out/$bcftools_out/${variant}_${bcftools_out}_normalized_dp${min_read_depth}.vcf',
    source = bcftools_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

bcftools_calls = env.Command(
    target = '$out/$checked_out/$bcftools_out/${variant}_${bcftools_out}_normalized_dp${min_read_depth}_checked.csv',
    source = [bcftools_cov_filtered,
              variant_cov_csv,
              bcftools_cov_csv],
    action = 'python $check_call $bcftools_out $variant $SOURCES $TARGET'
)


freebayes_cov_filtered = env.Command(
    target = '$out/$called_out/$freebayes_out/${variant}_${freebayes_out}_normalized_dp${min_read_depth}.vcf',
    source = freebayes_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

freebayes_calls = env.Command(
    target = '$out/$checked_out/$freebayes_out/${variant}_${freebayes_out}_normalized_dp${min_read_depth}_checked.csv',
    source = [freebayes_cov_filtered,
              variant_cov_csv,
              freebayes_cov_csv],
    action = 'python $check_call $freebayes_out $variant $SOURCES $TARGET'
)


deepvariant_cov_filtered = env.Command(
    target = '$out/$called_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_dp${min_read_depth}.vcf',
    source = deepvariant_pass,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

deepvariant_calls = env.Command(
    target = '$out/$checked_out/$deepvariant_out/${variant}_${deepvariant_out}_normalized_PASS_dp${min_read_depth}_checked.csv',
    source = [deepvariant_cov_filtered,
              variant_cov_csv,
              deepvariant_cov_csv],
    action = 'python $check_call $deepvariant_out $variant $SOURCES $TARGET'
)


discosnp_cov_filtered = env.Command(
    target = '$out/$called_out/$discosnp_out/${variant}_${discosnp_out}_normalized_dp${min_read_depth}.vcf',
    source = discosnp_pass_sorted,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

discosnp_calls = env.Command(
    target = '$out/$checked_out/$discosnp_out/${variant}_${discosnp_out}-edit_normalized_PASSsorted_dp${min_read_depth}_checked.csv',
    source = [discosnp_cov_filtered,
              variant_cov_csv,
              discosnp_cov_csv],
    action = 'python $check_call $discosnp_out $variant $SOURCES $TARGET'
)


lancet_cov_filtered = env.Command(
    target = '$out/$called_out/$lancet_out/${variant}_${lancet_out}_normalized_dp${min_read_depth}.vcf',
    source = lancet_final_vcf,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

lancet_calls = env.Command(
    target = '$out/$checked_out/$lancet_out/${variant}_${lancet_out}_normalized_dp${min_read_depth}_checked.csv',
    source = [lancet_cov_filtered,
              variant_cov_csv,
              lancet_cov_csv],
    action = 'python $check_call $lancet_out $variant $SOURCES $TARGET'
)


vardict_cov_filtered = env.Command(
    target = '$out/$called_out/$vardict_out/${variant}_${vardict_out}_normalized_dp${min_read_depth}.vcf',
    source = vardict_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' -o $TARGET $SOURCE'
)

vardict_calls = env.Command(
    target = '$out/$checked_out/$vardict_out/${variant}_${vardict_out}_normalized_dp${min_read_depth}_checked.csv',
    source = [vardict_cov_filtered,
              variant_cov_csv,
              vardict_cov_csv],
    action = 'python $check_call $vardict_out $variant $SOURCES $TARGET'
)

all_checked_csv = env.Command(
    target = '$out/$checked_out/${variant}_alltools_normalized_dp${min_read_depth}_checked.csv',
    source = [gatk_calls,
              bcftools_calls,
              freebayes_calls,
              deepvariant_calls,
              discosnp_calls,
              lancet_calls,
              vardict_calls],
    action = ('echo \'CHROM,POS,REF,ALT,CALLED_TYPE,TYPE,INS_TYPE,QUAL,AD_REF,AD_ALT,DP,BAM_DP,GT,ZYG,RK_DISCOSNP,TOOL,SAMPLE,TRUE_POS,FALSE_POS,FALSE_NEG\' > $TARGET; '
              'cat $SOURCES | sed \'/CHROM,POS,REF,ALT,CALLED_TYPE,TYPE,INS_TYPE,QUAL,AD_REF,AD_ALT,DP,BAM_DP,GT,ZYG,RK_DISCOSNP,TOOL,SAMPLE,TRUE_POS,FALSE_POS,FALSE_NEG/d\' >> $TARGET')
)


