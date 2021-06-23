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
    ENV = dict(os.environ, PATH=PATH, SHELLOPTS='errexit:pipefail'),
    variables = vars,
    SHELL = 'bash',
    cwd = os.getcwd(),
    venv_config = config.get('DEFAULT', 'virtualenv'),
    reference = config.get('DEFAULT', 'reference_genome'),
    ref_name = config.get('DEFAULT', 'reference_name'),
    accession = config.get('DEFAULT', 'reference_accession'),
    variant = config.get('variant_simulation', 'variant_name'),
    variants_config = config.get('variant_simulation', 'variants_config'),
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
    freebayes_out = config.get('variant_calling', 'freebayes_output'),
    discosnp_out = config.get('variant_calling', 'discosnp_output'),
    deepvariant_out = config.get('variant_calling', 'deepvariant_output'),
    vardict_scripts = config.get('variant_calling', 'vardict_scripts'),
    vardict_out = config.get('variant_calling', 'vardict_output'),
    delly_out = config.get('variant_calling', 'delly_output'),
    lancet_out = config.get('variant_calling', 'lancet_output'),
    ploidy = config.get('variant_calling', 'ploidy'),
    allele_fraction = config.get('variant_calling', 'min_allele_fraction'),
    min_read_depth = config.get('variant_calling', 'min_read_depth'),
    max_read_depth = config.get('variant_calling', 'max_read_depth'),
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
    gatk = '{} exec -B $cwd {} gatk'.format(singularity, gatk_img),
    cutadapt = '{} exec -B $cwd {} cutadapt'.format(singularity, cutadapt_img),
    bwa = '{} exec -B $cwd {} bwa mem'.format(singularity, bwa_img),
    samtools = '{} exec -B $cwd {} samtools'.format(singularity, samtools_img),
    bcftools = '{} exec -B $cwd {} bcftools'.format(singularity, bcftools_img),
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
    target = ['$out/$variants_out/${variant}.txt', '$out/$variants_out/${variant}.fa'],
    source = '$reference',
    action = ('python bin/variants.py --settings $variants_config $SOURCE $TARGETS')
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
              '> ${TARGETS[-1]} 2<&1')
)

# ########### Simulate NGS Reads ##############
sim_R1, sim_R2, simreads_log = env.Command(
    target = ['$out/$reads_out/${variant}_R1.fq', 
              '$out/$reads_out/${variant}_R2.fq', 
              '$log/$reads_out/${variant}_art_illumina.log'],
    source = simulated_variants_fa,
    action = ('./${venv_config}/bin/art_illumina -1 $read1_q -2 $read2_q -p -sam '
              '-i $SOURCE -l $read_len -f $read_depth -m $frag_len -s $sd '
              '-o $out/$reads_out/${variant}_R > ${TARGETS[-1]} 2<&1')
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
              '> ${TARGETS[-1]} 2<&1')
)

# ################ Map Reads ###################

sam = env.Command(
    target = ['$out/$mapped_out/${variant}_trimmed.sam'],
    source = ['$reference', R1trimmed, R2trimmed],
    action = ('$bwa $SOURCES -K $bwa_k '
              '-R \'@RG\\tID:${variant}\\tLB:LB_${variant}\\tPL:${rg_pl}\\tPU:${rg_pu}\\tSM:${variant}\' '
              '> $TARGET')
)

sorted_bam = env.Command(
    target = ['$out/$mapped_out/${variant}_trimmed-sorted.bam'],
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

# ################### Filter Reads ####################

mq_filtered_bam = env.Command(
    target = '$out/$deduped_out/${variant}_deduped_mq.bam',
    source = deduped,
    action = '$samtools view $SOURCE -q $mapq -bo $TARGET'
)

indexed_bam = env.Command(
    target = '$out/$deduped_out/${variant}_deduped_mq.bam.bai',
    source = mq_filtered_bam,
    action = '$samtools index $SOURCE'
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
              '$log/$called_out/${variant}_${gatk_out}_haplotypecaller.log'],
    source = ['$reference',
              mq_filtered_bam],
    action = ('$gatk HaplotypeCaller -ploidy $ploidy -R ${SOURCES[0]} -I ${SOURCES[1]} '
              '-O ${TARGETS[0]} -ERC GVCF > ${TARGETS[-1]} 2>&1')
)

gatk_gt, gatk_gt_log = env.Command(
    target = ['$out/$called_out/$gatk_out/${variant}_${gatk_out}.vcf',
              '$log/$called_out/${variant}_${gatk_out}_genotypegvcfs.log'],
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
              '$log/$called_out/${variant}_${gatk_out}_normalized.log'],
    source = ['$reference',
              gatk_gt],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

# ##################### bcftools #######################

pileup = env.Command(
    target = '$out/$called_out/$bcftools_out/${variant}_pileup.vcf',
    source = ['$reference',
              mq_filtered_bam],
    action = ('$samtools mpileup -m $min_read_depth -F $allele_fraction -u -f ${SOURCES[0]} '
              '-d $max_read_depth -A -B ${SOURCES[1]} -vo $TARGET')
)

bcftools_gvcf = env.Command(
    target = '$out/$called_out/$gvcf_out/${variant}_${bcftools_out}.g.vcf',
    source = pileup,
    action = '$bcftools call -g $min_read_depth -mO v -o $TARGET $SOURCE' 
)

bcftools_vcf = env.Command(
    target = '$out/$called_out/$bcftools_out/${variant}_${bcftools_out}.vcf',
    source = bcftools_gvcf,
    action = '$bcftools convert --gvcf2vcf $SOURCE -O v -o $TARGET'
)

bcftools_g_normalized, bcftools_g_norm_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${bcftools_out}_normalized.g.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${bcftools_out}_normalized.g.log'],
    source = ['$reference',
              bcftools_gvcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

bcftools_normalized, bcftools_norm_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${bcftools_out}_normalized.vcf',
              '$log/$called_out/$gvcf_out/${variant}_${bcftools_out}_normalized.log'],
    source = ['$reference',
              bcftools_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
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
    target = ['$out/$called_out/${variant}_${freebayes_out}_normalized.vcf',
              '$log/$called_out/${variant}_${freebayes_out}_normalized.log'],
    source = ['$reference',
              freebayes_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)


# ################# DeepVariant #################

deepvariant_gvcf, deepvariant_vcf, deepvariant_log = env.Command(
    target = ['$out/$called_out/$gvcf_out/${variant}_${deepvariant_out}.g.vcf',
              '$out/$called_out/${variant}_${deepvariant_out}.vcf',
              '$log/$called_out/${variant}_${deepvariant_out}.log'],
    source = ['$reference',
              mq_filtered_bam],
    action = ('$deepvariant --model_type=WGS --ref=/input/${SOURCES[0]} --reads=/input/${SOURCES[1]} '
              '--output_gvcf=/output/${TARGETS[0]} --output_vcf=/output/${TARGETS[1]} --num_shards=$max_threads '
              '--logging_dir=/output/$out/$called_out/${deepvariant_out}/ > ${TARGETS[-1]} 2>&1')
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
    target = ['$out/$called_out/${variant}_${deepvariant_out}_normalized.vcf',
              '$log/$called_out/${variant}_${deepvariant_out}_normalized.log'],
    source = ['$reference',
              deepvariant_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

deepvariant_pass = env.Command(
    target = '$out/$called_out/${variant}_${deepvariant_out}_normalized_PASS.vcf',
    source = deepvariant_normalized,
    action = 'grep "#" $SOURCE > $TARGET; grep "$$(printf \'\\t\')PASS$$(printf \'\\t\')" $SOURCE >> $TARGET'
)

# ################### delly #####################

delly_bcf, delly_log  = env.Command(
    target = ['$out/$called_out/${variant}_${delly_out}.bcf',
              '$log/$called_out/${variant}_${delly_out}.log'],
    source = ['$reference',
              mq_filtered_bam],
    action = '$delly call -g ${SOURCES[0]} -o ${TARGETS[0]} ${SOURCES[1]} > ${TARGETS[-1]} 2>&1'
)

delly_vcf = env.Command(
    target = '$out/$called_out/${variant}_${delly_out}.vcf',
    source = delly_bcf,
    action = '$bcftools view $SOURCE -O v -o $TARGET'
)

delly_normalized, delly_norm_log = env.Command(
    target = ['$out/$called_out/${variant}_${delly_out}_normalized.vcf',
              '$log/$called_out/${variant}_${delly_out}_normalized.log'],
    source = ['$reference',
              delly_vcf],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
)

# ################### Lancet ####################

#lancet_vcf = env.Command(
#    target = '$out/$called_out/${variant}_${lancet_out}/.vcf',
#    source = ['$reference',
#              mq_filtered_bam,
#              '$out/$deduped_out/${ref_name}_deduped_mq10.bam'],
#    action = './bin/lancet --tumor ${SOURCES[1]} --normal ${SOURCES[2]} --ref ${SOURCES[0]} '
#             ' --reg $accession 
#)

# ################## DiscoSnp ###################

ref_fof = env.Command(
    target = '$out/$called_out/$discosnp_out/${ref_name}_fof.txt',
    source = None,
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

#discosnp_vcf, discosnp_log = env.Command(
#    target = ['$out/${ref_name}_${variant}_k_${kmer_size}_c_${coverage}_D_100_P_${snp_per_bubble}_b_${disco_mode}_coherent.vcf',
#              '$log/$called_out/${ref_name}_${variant}_discosnp.log'],
#    source = ['$reference',
#              '$out/$called_out/$discosnp_out/${ref_name}_${variant}_fof.txt'],
#    action = ('$discosnp $out -r $called_out/$discosnp_out/${ref_name}_${variant}_fof.txt -P $snp_per_bubble '
#              '-b $disco_mode -k $kmer_size -c $coverage -T -l '
#              '-G ../$reference -p ${ref_name}_${variant} -u $max_threads > ${TARGETS[-1]} 2>&1')
#)

#discosnp_normalized, discosnp_norm_log = env.Command(
#    target = ['$out/$called_out/$discosnp_out/${ref_name}_${variant}_discosnp_normalized.vcf',
#              '$log/$called_out/$discosnp_out/${ref_name}_${variant}_discosnp_normalized.log'],
#    source = ['$reference',
#              discosnp_vcf],
#    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
#              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
#)


# ################### VarDict ####################

bedfile = env.Command(
    target = '$out/$deduped_out/${variant}_deduped_mq10.bed',
    source = mq_filtered_bam,
    action = '$bedtools bamtobed -i $SOURCE > $TARGET'
)

#vardict_vcf = env.Command(
#    target = '$out/$called_out/${variant}_vardict.vcf',
#    source = ['$reference',
#              mq_filtered_bam,
#              bedfile],
#    action = ('$vardict -G ${SOURCES[0]} -f $allele_fraction -N $variant -b ${SOURCES[1]} '
#              '-c 1 -S 2 -E 3 -g 4 ${SOURCES[2]} | ${vardict_scripts}/teststrandbias.R '
#              '| ${vardict_scripts}/var2vcf_valid.pl -N $variant -E -f $allele_fraction '
#              '> $TARGET')
#)

#vardict_normalized, vardict_norm_log = env.Command(
#    target = ['$out/$called_out/${variant}_${vardict_out}_normalized.vcf',
#              '$log/$called_out/${variant}_${vardict_out}_normalized.log'],
#    source = ['$reference',
#              vardict_vcf],
#    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]} '
#              '-O ${TARGETS[0]} > ${TARGETS[-1]} 2>&1')
#)
