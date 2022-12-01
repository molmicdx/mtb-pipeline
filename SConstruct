import os 
import sys
from bin.utils import parse_config

SETTINGS = ARGUMENTS.get('SETTINGS', 'configs/settings.conf')
config = parse_config(SETTINGS)

# Ensure that a virtualenv is active before importing non-stdlib dependencies.
try:
    venv = os.environ.get('VIRTUAL_ENV')
except KeyError:
    sys.exit('an active virtualenv is required')

singularity = config.get('singularity', 'singularity')
art_img = config.get('singularity', 'art')
gatk_img = config.get('singularity', 'gatk')
cutadapt_img = config.get('singularity', 'cutadapt')
seqmagick_img = config.get('singularity', 'seqmagick')
bwa_img = config.get('singularity', 'bwa')
samtools_img = config.get('singularity', 'samtools')
bcftools_img = config.get('singularity', 'bcftools')
bedtools_img = config.get('singularity', 'bedtools')
discosnp_img = config.get('singularity', 'discosnp')
freebayes_img = config.get('singularity', 'freebayes')
deepvariant_img = config.get('singularity', 'deepvariant')
vardict_img = config.get('singularity', 'vardict')
lancet_img = config.get('singularity', 'lancet')
htslib_img = config.get('singularity', 'htslib')
igv_img = config.get('singularity', 'igv-reports')
py_img = config.get('singularity', 'py-deps')

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
    ENV = dict(os.environ, PATH=PATH),
    variables = vars,
    SHELL = 'bash',
    cwd = os.getcwd(),
    venv_config = config.get('DEFAULT', 'virtualenv'),
    ref_dir = config.get('DEFAULT', 'reference_dir'),
    reference = config.get('DEFAULT', 'reference_genome'),
    ref_name = config.get('DEFAULT', 'reference_name'),
    accession = config.get('DEFAULT', 'reference_accession'),
    sample = config.get('variant_simulation', 'variant_name'),
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
    gatk_out = config.get('variant_calling', 'gatk_output'),
    bcftools_out = config.get('variant_calling', 'bcftools_output'),
    ploidy_file = config.get('variant_calling', 'bcftools_ploidy_file'),
    bcftools_ann = config.get('variant_calling', 'bcftools_annotation'),
    freebayes_out = config.get('variant_calling', 'freebayes_output'),
    discosnp_out = config.get('variant_calling', 'discosnp_output'),
    deepvariant_out = config.get('variant_calling', 'deepvariant_output'),
    vardict_scripts = config.get('variant_calling', 'vardict_scripts'),
    vardict_out = config.get('variant_calling', 'vardict_output'),
    lancet_out = config.get('variant_calling', 'lancet_output'),
    min_reads = config.get('variant_calling', 'min_reads'),
    max_reads = config.get('variant_calling', 'max_reads'),
    rg_pl = config.get('read_mapping', 'read_group_PL'),
    rg_pu = config.get('read_mapping', 'read_group_PU'),
    mapq = config.get('read_mapping', 'min_mapq'),
    min_read_depth = config.get('report', 'filter_read_depth'),
    mutation_to_flank = config.get('report', 'ignore_cov_on_mutation'), 
    num_flanking_bp = config.get('report', 'num_flanking_bp_to_use'),
    igv_info = config.get('report', 'igv_info'),
    igv_flank = config.get('report', 'igv_flanking_sites'),
    art = '{} exec -B $cwd {} art_illumina'.format(singularity, art_img),
    gatk = '{} exec -B $cwd {} gatk'.format(singularity, gatk_img),
    cutadapt = '{} exec -B $cwd {} cutadapt'.format(singularity, cutadapt_img),
    seqmagick = '{} run -B $cwd {} seqmagick'.format(singularity, seqmagick_img),
    bwa = '{} exec -B $cwd {} bwa'.format(singularity, bwa_img),
    samtools = '{} exec -B $cwd {} samtools'.format(singularity, samtools_img),
    bcftools = '{} exec -B $cwd {}'.format(singularity, bcftools_img),
    bedtools = '{} exec -B $cwd {} bedtools'.format(singularity, bedtools_img),
    discosnp = '{} run --pwd $cwd -B $cwd {}'.format(singularity, discosnp_img),
    freebayes = '{} exec -B $cwd {} freebayes'.format(singularity, freebayes_img),
    deepvariant = '{} exec -B $cwd {} python3 /opt/deepvariant/bin/run_deepvariant.py' \
                  .format(singularity, deepvariant_img),
    lancet = '{} exec -B $cwd {} lancet'.format(singularity, lancet_img),
    vardict = '{} exec -B $cwd {} vardict-java'.format(singularity, vardict_img),
    bgzip = '{} exec -B $cwd {} bgzip'.format(singularity, htslib_img),
    tabix = '{} exec -B $cwd {} tabix'.format(singularity, htslib_img),
    igv = '{} exec -B $cwd {} create_report'.format(singularity, igv_img),
    py_deps = '{} exec -B $cwd {} python3'.format(singularity, py_img)
)


# Help(vars.GenerateHelpText(env))
# ############### start inputs ################

# ######### Index Reference Sequence #########

fa_index = env.Command(
    target = '${ref_dir}/${reference}.fa.fai',
    source = '${ref_dir}/${reference}.fa',
    action = '$samtools faidx $SOURCE'
)

bwa_idx = env.Command(
    target = '${ref_dir}/${reference}.fa.bwt',
    source = '${ref_dir}/${reference}.fa',
    action = '$bwa index $SOURCE'
)

# ######## Create Reference Seq Dictionary ########

ref_dict = env.Command(
    target = '${ref_dir}/${reference}.dict',
    source = '${ref_dir}/${reference}.fa',
    action = '$gatk CreateSequenceDictionary -R $SOURCE'
)

# ############# Simulate Variants #############
simulated_variants_table, simulated_variants_fa, svlog = env.Command(
    target = ['$out/${sample}/variants_table.txt',
              '$out/${sample}/variant_genome.fa',
              '$log/${sample}/variants.log'], 
    source = '${ref_dir}/${reference}.fa',
    action = ('$py_deps bin/variants.py --settings $variants_config $SOURCE '
              '${TARGETS[0]} ${TARGETS[1]} > ${TARGETS[-1]} 2>&1 ')
)

# ############# Normalize Variant VCF ###############
simulated_variant_vcf = env.Command(
    target = '$out/${sample}/variants_table.txt.vcf',
    source = simulated_variants_table,
    action = '$py_deps bin/to_vcf.py $SOURCE $sample'
)

normalized_variant_vcf = env.Command(
    target = '$out/${sample}/variants_normalized.vcf', 
    source = simulated_variant_vcf,
    action = ('$gatk LeftAlignAndTrimVariants -R ${ref_dir}/${reference}.fa '
              '-V $SOURCE -O $TARGET > $log/${sample}/normalized.log 2>&1')
)

########### Simulate NGS Reads ##############

ref_R1, ref_R2 = env.Command(
    target = ['$out/${ref_name}/R1.fq',
              '$out/${ref_name}/R2.fq'],
    source = '${ref_dir}/${reference}.fa',
    action = ('$art -1 ${ref_dir}/${read1_q} -2 ${ref_dir}/${read2_q} -p -sam '
              '-i $SOURCE -l $read_len -f $read_depth -m $frag_len -s $sd '
              '-o $out/${ref_name}/R')
)

gzrefR1, gzrefR2 = env.Command(
    target = ['$out/${ref_name}/R1.fq.gz',
              '$out/${ref_name}/R2.fq.gz'],
    source = [ref_R1, ref_R2],
    action = ('gzip $SOURCES --keep')
)

sim_R1, sim_R2 = env.Command(
    target = ['$out/${sample}/R1.fq', 
              '$out/${sample}/R2.fq'],
    source = simulated_variants_fa,
    action = ('$art -1 ${ref_dir}/${read1_q} -2 ${ref_dir}/${read2_q} -p -sam '
              '-i $SOURCE -l $read_len -f $read_depth -m $frag_len -s $sd '
              '-o $out/${sample}/R > $log/${sample}/art_illumina.log 2>&1')
)

gzsimR1, gzsimR2 = env.Command(
    target = ['$out/${sample}/R1.fq.gz',
              '$out/${sample}/R2.fq.gz'],
    source = [sim_R1, sim_R2],
    action = ('gzip $SOURCES --keep')
)


# ############## Trim NGS Reads ################
refR1trimmed, refR2trimmed = env.Command(
    target = ['$out/${ref_name}/R1.trimmed.fq.gz',
              '$out/${ref_name}/R2.trimmed.fq.gz'],
    source = [gzrefR1, gzrefR2],
    action = ('$cutadapt --cores 8 -q $min_read_q -b file:${ref_dir}/${adaptors} -B '
              'file:${ref_dir}/${adaptors} --minimum-length $min_read_len '
              '-o ${TARGETS[0]} -p ${TARGETS[1]} $SOURCES ')
)

R1trimmed, R2trimmed = env.Command(
    target = ['$out/${sample}/R1.trimmed.fq.gz',
              '$out/${sample}/R2.trimmed.fq.gz'],
    source = [gzsimR1, gzsimR2],
    action = ('$cutadapt --cores 8 -q $min_read_q -b file:${ref_dir}/${adaptors} -B '
              'file:${ref_dir}/${adaptors} --minimum-length $min_read_len '
              '-o ${TARGETS[0]} -p ${TARGETS[1]} $SOURCES '
              '> $log/${sample}/cutadapt.log 2>&1')
)

# ############ Get Reads Info #################

seq_log = env.Command(
    target = '$log/${sample}/seqmagick.log',
    source = [gzsimR1, gzsimR2, 
              R1trimmed, R2trimmed],
    action = '$seqmagick info $SOURCES > $TARGET'
) 


# ################ Map Reads ###################

ref_sam = env.Command(
    target = '$out/${ref_name}/trimmed.sam',
    source = ['${ref_dir}/${reference}.fa',
              refR1trimmed, refR2trimmed],
    action = ('$bwa mem $SOURCES -K 100000000 -R '
	      '\'@RG\\tID:${accession}\\tLB:LB_${accession}\\tPL:${rg_pl} \\tPU:${rg_pu}\\tSM:${accession}\' > $TARGET ')
)

ref_sorted_bam = env.Command(
    target = '$out/${ref_name}/trimmed-sorted.bam',
    source = ref_sam,
    action = '$samtools sort $SOURCE $out/${ref_name}/trimmed-sorted'
)

sam = env.Command(
    target = '$out/${sample}/trimmed.sam',
    source = ['${ref_dir}/${reference}.fa',
              R1trimmed, R2trimmed],
    action = ('$bwa mem $SOURCES -K 100000000 -R '
	      '\'@RG\\tID:${sample}\\tLB:LB_${sample}\\tPL:${rg_pl} \\tPU:${rg_pu}\\tSM:${sample}\' > $TARGET '
              '2> $log/${sample}/trimmed_mapped.log')
)

sorted_bam = env.Command(
    target = '$out/${sample}/trimmed-sorted.bam',
    source = sam,
    action = '$samtools sort $SOURCE $out/${sample}/trimmed-sorted'
)

# ############## Remove Duplicate Reads ###############
ref_deduped = env.Command(
    target = '$out/${ref_name}/deduped.bam', 
    source = ref_sorted_bam,
    action = ('$gatk MarkDuplicates -I $SOURCE -O $TARGET '
              '-M $out/${ref_name}/deduped_metrics.txt '
              '--REMOVE_DUPLICATES TRUE')
)

deduped = env.Command(
    target = '$out/${sample}/deduped.bam', 
    source = sorted_bam,
    action = ('$gatk MarkDuplicates -I $SOURCE -O $TARGET '
              '-M $out/${sample}/deduped_metrics.txt '
              '--REMOVE_DUPLICATES TRUE > $log/${sample}/deduped.log 2>&1')
)

# ################### Filter Reads, Get BAM Read Depths ####################

ref_mq_filtered_bam = env.Command(
    target = '$out/${ref_name}/deduped_mq.bam',
    source = ref_deduped,
    action = ('$samtools view $SOURCE -q $mapq -bo $TARGET; '
              '$samtools index $TARGET')
)

ref_genome_cov = env.Command(
    target = '$out/${ref_name}/deduped_mq_genomecov.bed',
    source = ref_mq_filtered_bam,
    action = '$bedtools genomecov -ibam $SOURCE -bga > $TARGET'
)
mq_filtered_bam = env.Command(
    target = '$out/${sample}/deduped_mq.bam',
    source = deduped,
    action = ('$samtools view $SOURCE -q $mapq -bo $TARGET; '
              '$samtools index $TARGET')
)

genome_cov = env.Command(
    target = '$out/${sample}/deduped_mq_genomecov.bed',
    source = mq_filtered_bam,
    action = '$bedtools genomecov -ibam $SOURCE -bga > $TARGET'
)

meandepth = env.Command(
    target = ['$log/${sample}/meandepth.log'],
    source = [sorted_bam, mq_filtered_bam],
    action = ('./bin/get_coverage.sh ${SOURCES[0]} $accession > $TARGET; '
              './bin/get_coverage.sh ${SOURCES[1]} $accession >> $TARGET')
)


# ############### Get Variant Coverage ################

variant_formatted_csv, variant_bed = env.Command(
    target = ['$out/${sample}/formatted.csv',
              '$out/${sample}/variant.bed'],
    source = normalized_variant_vcf,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

variant_cov_bed, variant_cov_csv = env.Command(
    target = ['$out/${sample}/genomecov_intersect.bed',
              '$out/${sample}/variant_cov.csv'],
    source = [variant_formatted_csv,
              variant_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} '
              '-wo > ${TARGETS[0]}; '
	      '$py_deps bin/add_cov.py ${TARGETS[0]} ${SOURCES[0]} '
              '${TARGETS[1]} --sample $sample')
)


# ############### end inputs ##################

# ################# Call Variants #####################

# ##################### GATK ##########################

gatk_gvcf = env.Command(
    target = '$out/${sample}/${gatk_out}/sample.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$gatk HaplotypeCaller -ploidy 1 -R ${SOURCES[0]} '
              '-I ${SOURCES[1]} -O $TARGET -ERC GVCF')
)

gatk_gt = env.Command(
    target = '$out/${sample}/${gatk_out}/sample.vcf',
    source = ['${ref_dir}/${reference}.fa',
              gatk_gvcf],
    action = ('$gatk GenotypeGVCFs -R ${SOURCES[0]} -V ${SOURCES[1]} '
              '-O $TARGET')
)

gatk_g_normalized = env.Command(
    target = '$out/${sample}/${gatk_out}/sample_normalized.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              gatk_gvcf, 
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

gatk_normalized = env.Command(
    target = '$out/${sample}/${gatk_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              gatk_gt,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

gatk_bgz = env.Command(
    target = '$out/${sample}/${gatk_out}/sample_normalized.vcf.gz',
    source = gatk_normalized,
    action = '$bgzip < $SOURCE > $TARGET'
)

gatk_tbi, gatk_igv = env.Command(
    target = ['$out/${sample}/${gatk_out}/sample_normalized.vcf.gz.tbi',
              '$out/${sample}/${gatk_out}/igv.html'],
    source = [gatk_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; $igv ${SOURCES[0]} ${SOURCES[1]} '
              '--flanking $igv_flank --info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

gatk_csv, gatk_bed = env.Command(
    target = ['$out/${sample}/${gatk_out}/sample_normalized.vcf.csv',
              '$out/${sample}/${gatk_out}/sample_normalized.vcf.csv.bed'],
    source = gatk_normalized,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank ' 
              '--bp $num_flanking_bp')
)

gatk_cov_bed, gatk_cov_csv = env.Command(
    target = ['$out/${sample}/${gatk_out}/'
              'sample_normalized_genomecov_intersect.bed',
              '$out/${sample}/${gatk_out}/sample_normalized_cov.csv'],
    source = [gatk_csv,
              gatk_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; $py_deps bin/add_cov.py ${TARGETS[0]} '
              '${SOURCES[0]} ${TARGETS[1]}')
)


# ##################### bcftools #######################

bcftools_gvcf = env.Command(
    target = '$out/${sample}/${bcftools_out}/sample.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$bcftools bcftools mpileup -Ou -f $SOURCES -d $max_reads '
              '-a $bcftools_ann | '
	      '$bcftools bcftools call -Ov -m --gvcf 1 --ploidy-file '
              '${ref_dir}/${ploidy_file} -o $TARGET')
)

bcftools_g_normalized = env.Command(
    target = '$out/${sample}/${bcftools_out}/sample_normalized.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              bcftools_gvcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

bcftools_vcf = env.Command(
    target = '$out/${sample}/${bcftools_out}/sample.vcf',
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$bcftools bcftools mpileup -Ou -f $SOURCES -d $max_reads '
              '-a $bcftools_ann | '
	      '$bcftools bcftools call -Ov -mv --ploidy-file ${ref_dir}/${ploidy_file} '
              '-o $TARGET')
)

bcftools_normalized = env.Command(
    target = '$out/${sample}/${bcftools_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              bcftools_vcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

bcftools_bgz = env.Command(
    target = '$out/${sample}/${bcftools_out}/sample_normalized.vcf.gz',
    source = bcftools_normalized,
    action = '$bgzip < $SOURCE > $TARGET'
)

bcftools_tbi, bcftools_igv = env.Command(
    target = ['$out/${sample}/${bcftools_out}/sample_normalized.vcf.gz.tbi',
              '$out/${sample}/${bcftools_out}/igv.html'],
    source = [bcftools_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; '
	      '$igv ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank '
              '--info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

bcftools_csv, bcftools_bed = env.Command(
    target = ['$out/${sample}/${bcftools_out}/sample_normalized.vcf.csv',
              '$out/${sample}/${bcftools_out}/sample_normalized.vcf.csv.bed'],
    source = bcftools_normalized,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

bcftools_cov_bed, bcftools_cov_csv = env.Command(
    target = ['$out/${sample}/${bcftools_out}/'
              'sample_normalized_genomecov_intersect.bed',
              '$out/${sample}/${bcftools_out}/sample_normalized_cov.csv'],
    source = [bcftools_csv,
              bcftools_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; $py_deps bin/add_cov.py ${TARGETS[0]} '
              '${SOURCES[0]} ${TARGETS[1]}')
)


# ##################### FreeBayes  #######################

freebayes_gvcf = env.Command(
    target = '$out/${sample}/${freebayes_out}/sample.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$freebayes -f ${SOURCES[0]} -p 1 --min-alternate-fraction 0.2 '
              '--gvcf ${SOURCES[1]} > $TARGET')
)

fb_allgt_vcf, freebayes_vcf = env.Command(
    target = ['$out/${sample}/${freebayes_out}/sample_allgt.vcf',
              '$out/${sample}/${freebayes_out}/sample.vcf'],
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$freebayes -f ${SOURCES[0]} -p 1 --min-alternate-fraction 0.2 '
              '${SOURCES[1]} > ${TARGETS[0]}; '
	      '$bcftools bcftools plugin setGT ${TARGETS[0]} -- -t q -n . '
              '-i \'GT="0"\' > tmp.vcf; '
              '$bcftools bcftools view -g ^miss tmp.vcf > ${TARGETS[1]}')
)

freebayes_g_normalized = env.Command(
    target = '$out/${sample}/${freebayes_out}/sample_normalized.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              freebayes_gvcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

freebayes_normalized = env.Command(
    target = '$out/${sample}/${freebayes_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              freebayes_vcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

freebayes_bgz = env.Command(
    target = '$out/${sample}/${freebayes_out}/sample_normalized.vcf.gz',
    source = freebayes_normalized,
    action = '$bgzip < $SOURCE > $TARGET'
)

freebayes_tbi, freebayes_igv = env.Command(
    target = ['$out/${sample}/${freebayes_out}/sample_normalized.vcf.gz.tbi',
              '$out/${sample}/${freebayes_out}/igv.html'],
    source = [freebayes_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; '
	      '$igv ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank '
              '--info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

freebayes_csv, freebayes_bed = env.Command(
    target = ['$out/${sample}/${freebayes_out}/sample_normalized.vcf.csv',
              '$out/${sample}/${freebayes_out}/sample_normalized.vcf.csv.bed'],
    source = freebayes_normalized,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

freebayes_cov_bed, freebayes_cov_csv = env.Command(
    target = ['$out/${sample}/${freebayes_out}/'
              'sample_normalized_genomecov_intersect.bed',
              '$out/${sample}/${freebayes_out}/sample_normalized_cov.csv'],
    source = [freebayes_csv,
              freebayes_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; '
              '$py_deps bin/add_cov.py ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################# DeepVariant #################

deepvariant_gvcf, deepvariant_vcf = env.Command(
    target = ['$out/${sample}/${deepvariant_out}/sample.g.vcf',
              '$out/${sample}/${deepvariant_out}/sample.vcf'],
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$deepvariant --model_type=WGS --ref=${SOURCES[0]} '
              '--reads=${SOURCES[1]} --output_gvcf=${TARGETS[0]} '
              '--output_vcf=${TARGETS[1]} --num_shards=12 '
              '--logging_dir=$out/${sample}/${deepvariant_out}/$log')
)

deepvariant_g_normalized = env.Command(
    target = '$out/${sample}/${deepvariant_out}/sample_normalized.g.vcf',
    source = ['${ref_dir}/${reference}.fa',
              deepvariant_gvcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

deepvariant_normalized = env.Command(
    target = '$out/${sample}/${deepvariant_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              deepvariant_vcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

deepvariant_pass = env.Command(
    target = '$out/${sample}/${deepvariant_out}/sample_normalized_PASS.vcf',
    source = deepvariant_normalized,
    action = ('grep "#" $SOURCE > $TARGET;'
              'grep "$$(printf \'\\t\')PASS$$(printf\'\\t\')" '
              '$SOURCE >> $TARGET')
)

deepvariant_bgz = env.Command(
    target = '$out/${sample}/${deepvariant_out}/sample_normalized_PASS.vcf.gz',
    source = deepvariant_pass,
    action = '$bgzip < $SOURCE > $TARGET'
)

deepvariant_tbi, deepvariant_igv = env.Command(
    target = ['$out/${sample}/${deepvariant_out}/'
              'sample_normalized_PASS.vcf.gz.tbi',
              '$out/${sample}/${deepvariant_out}/igv.html'],
    source = [deepvariant_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; '
	      '$igv ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank '
              '--info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

deepvariant_csv, deepvariant_bed = env.Command(
    target = ['$out/${sample}/${deepvariant_out}/sample_normalized.vcf.csv',
	      '$out/${sample}/${deepvariant_out}/'
              'sample_normalized.vcf.csv.bed'],
    source = deepvariant_pass,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

deepvariant_cov_bed, deepvariant_cov_csv = env.Command(
    target = ['$out/${sample}/${deepvariant_out}/'
              'sample_normalized_genomecov_intersect.bed',
              '$out/${sample}/${deepvariant_out}/sample_normalized_cov.csv'],
    source = [deepvariant_csv,
              deepvariant_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; '
	      '$py_deps bin/add_cov.py ${TARGETS[0]} ${SOURCES[0]} '
              '${TARGETS[1]}')
)


# ################### Lancet ####################

lancet_vcf = env.Command(
    target = '$out/${sample}/${lancet_out}/sample.vcf', 
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam,
              ref_mq_filtered_bam],
    action = ('$lancet --tumor ${SOURCES[1]} --normal ${SOURCES[2]} '
              '--ref ${SOURCES[0]} --reg ${accession} --min-vaf-tumor 0.2 '
              '--low-cov $min_read_depth --num-threads 12 > $TARGET '
              '2>/dev/null')
)

lancet_normalized = env.Command(
    target = '$out/${sample}/${lancet_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              lancet_vcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

lancet_pass = env.Command(
    target = '$out/${sample}/${lancet_out}/sample_normalized_PASS.vcf',
    source = lancet_normalized,
    action = ('grep "#" $SOURCE > $TARGET; '
              'grep "$$(printf \'\\t\')PASS$$(printf\'\\t\')" '
              '$SOURCE >> $TARGET')
)

lancet_final_vcf = env.Command(
    target = '$out/${sample}/${lancet_out}/sample_normalized_PASSsorted.vcf',
    source = lancet_pass,
    action = ('for sample in $$(zgrep -m 1 "^#CHROM" $SOURCE | cut -f10-); do '
	      '    $bcftools bcftools view -c 1 -Ov -s $$sample '
              '-o $out/$$sample\'_normalized_PASS.vcf\' $SOURCE; done; '
	      '$gatk SortVcf -I $out/${sample}_normalized_PASS.vcf -O $TARGET;'
              'rm $out/*_normalized_PASS.vcf')
)

lancet_bgz = env.Command(
    target = '$out/${sample}/${lancet_out}/sample_normalized_PASSsorted.vcf.gz',
    source = lancet_final_vcf,
    action = '$bgzip < $SOURCE > $TARGET'
)

lancet_tbi, lancet_igv = env.Command(
    target = ['$out/${sample}/${lancet_out}/'
              'sample_normalized_PASSsorted.vcf.gz.tbi',
              '$out/${sample}/${lancet_out}/igv.html'],
    source = [lancet_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; '
	      '$igv ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank '
              '--info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

lancet_csv, lancet_bed = env.Command(
    target = ['$out/${sample}/${lancet_out}/'
              'sample_normalized_PASSsorted.vcf.csv',
	      '$out/${sample}/${lancet_out}/'
              'sample_normalized_PASSsorted.vcf.csv.bed'],
    source = lancet_final_vcf,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

lancet_cov_bed, lancet_cov_csv = env.Command(
    target = ['$out/${sample}/${lancet_out}/'
              'sample_normalized_PASSsorted_genomecov_intersect.bed',
	      '$out/${sample}/${lancet_out}/'
              'sample_normalized_PASSsorted_cov.csv'],
    source = [lancet_csv,
              lancet_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; '
	      '$py_deps bin/add_cov.py ${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}')
)


# ################## DiscoSnp ###################

ref_fof = env.Command(
    target = '$out/${ref_name}/fof.txt',
    source = [refR1trimmed, 
              refR2trimmed],
    action = ('echo ${cwd}/${out}/${ref_name}/R1.trimmed.fq.gz > $TARGET; '
              'echo ${cwd}/${out}/${ref_name}/R2.trimmed.fq.gz >> $TARGET')
)

variant_fof = env.Command(
    target = '$out/${sample}/${discosnp_out}/sample_fof.txt',
    source = None,
    action = ('echo ${cwd}/${out}/${sample}/R1.trimmed.fq.gz > $TARGET; '
              'echo ${cwd}/${out}/${sample}/R2.trimmed.fq.gz >> $TARGET')
)

fof = env.Command(
    target = '$out/${sample}/${discosnp_out}/fof.txt',
    source = [ref_fof,
              variant_fof],
    action = ('echo ${cwd}/${out}/${ref_name}/fof.txt > $TARGET; '
              'echo sample_fof.txt >> $TARGET')
)

discosnp_vcf = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_k_31_c_auto_D_100_P_6_b_1_coherent.vcf',
    source = ['${ref_dir}/${reference}.fa',
              fof],
    action = ('$discosnp $out -r ../${SOURCES[1]} -P 6 '
              '-b 1 -k 31 -c auto -T -l '
              '-G ../${SOURCES[0]} -p sample -u 12; '
              'mv $out/sample* $out/$sample/$discosnp_out/; '
	      'sed -i \'s/INDEL_.*_path_[0-9]*/${accession}/g\' $TARGET') 
    #temp fix for inexplicable VCF entry; may falsely inflate false pos calls
)

discosnp_normalized = env.Command(
    target = '$out/${sample}/${discosnp_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              discosnp_vcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} -V ${SOURCES[1]}'
              ' -O $TARGET')
)

discosnp_formatted = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'ref_sample_discosnp-edit_normalized.vcf',
    source = discosnp_normalized,
    action = ('sed -e \'s/##SAMPLE/##sample/\' -e \'s/G1/${accession}/\' '
              '-e \'s/G2/${sample}/\' < $SOURCE > $TARGET')
)

discosnp_final_vcf = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_discosnp-edit_normalized.vcf',
    source = discosnp_formatted,
    action = ('for sample in $$(zgrep -m 1 "^#CHROM" $SOURCE | cut -f10-); do '
	      '    $bcftools bcftools view -c 1 -Ov -s $$sample '
              '-o $out/$$sample\'_discosnp-edit_normalized.vcf\' $SOURCE; done;' 
              'mv $out/${sample}_discosnp-edit_normalized.vcf $TARGET; ' 
              'rm $out/*_discosnp-edit_normalized.vcf')
)

discosnp_pass = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_discosnp-edit_normalized_PASS.vcf',
    source = discosnp_final_vcf,
    action = ('grep "#" $SOURCE > $TARGET; '
              'grep "$$(printf \'\\t\')PASS$$(printf\'\\t\')" '
              '$SOURCE >> $TARGET')
)

discosnp_pass_sorted = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_discosnp-edit_normalized_PASSsorted.vcf',
    source = discosnp_pass,
    action = '$gatk SortVcf -I $SOURCE -O $TARGET'
)

discosnp_bgz = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_discosnp-edit_normalized_PASSsorted.vcf.gz',
    source = discosnp_pass_sorted,
    action = '$bgzip < $SOURCE > $TARGET'
)

discosnp_tbi, discosnp_igv = env.Command(
    target = ['$out/${sample}/${discosnp_out}/'
              'sample_discosnp-edit_normalized_PASSsorted.vcf.gz.tbi',
              '$out/${sample}/${discosnp_out}/igv.html'],
    source = [discosnp_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; '
	      '$igv ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank '
              '--info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

discosnp_csv, discosnp_bed = env.Command(
    target = ['$out/${sample}/${discosnp_out}/'
              'sample_discosnp-edit_normalized_PASSsorted.vcf.csv',
	      '$out/${sample}/${discosnp_out}/'
              'sample_discosnp-edit_normalized_PASSsorted.vcf.csv.bed'],
    source = discosnp_pass_sorted,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

discosnp_cov_bed, discosnp_cov_csv = env.Command(
    target = ['$out/${sample}/${discosnp_out}/'
              'sample_discosnp-edit_normalized_PASSsorted_genomecov_intersect.bed',
	      '$out/${sample}/${discosnp_out}/'
              'sample_discosnp-edit_normalized_PASSsorted_cov.csv'],
    source = [discosnp_csv,
              discosnp_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; '
	      '$py_deps bin/add_cov.py ${TARGETS[0]} ${SOURCES[0]} '
              '${TARGETS[1]}')
)


# ################### VarDict ####################

bedfile = env.Command(
    target = '$out/${sample}/${vardict_out}/deduped_mq.bed',
    source = mq_filtered_bam,
    action = '$bedtools bamtobed -i $SOURCE > $TARGET'
)

vardict_vcf = env.Command(
    target = '$out/${sample}/${vardict_out}/sample.vcf',
    source = ['${ref_dir}/${reference}.fa',
              mq_filtered_bam,
              bedfile],
    action = ('$vardict -G ${SOURCES[0]} -f 0.2 -N $sample -b ${SOURCES[1]} '
	      '-c 1 -S 2 -E 3 -g 4 ${SOURCES[2]} | '
              '${vardict_scripts}/teststrandbias.R '
              '| ${vardict_scripts}/var2vcf_valid.pl -N $sample '
              '-E -f 0.2 > $TARGET')
)

vardict_normalized = env.Command(
    target = '$out/${sample}/${vardict_out}/sample_normalized.vcf',
    source = ['${ref_dir}/${reference}.fa',
              vardict_vcf,
              ref_dict],
    action = ('$gatk LeftAlignAndTrimVariants -R ${SOURCES[0]} '
              '-V ${SOURCES[1]} -O $TARGET')
)

vardict_bgz = env.Command(
    target = '$out/${sample}/${vardict_out}/sample_normalized.vcf.gz',
    source = vardict_normalized,
    action = '$bgzip < $SOURCE > $TARGET'
)

vardict_tbi, vardict_igv = env.Command(
    target = ['$out/${sample}/${vardict_out}/sample_normalized.vcf.gz.tbi',
              '$out/${sample}/${vardict_out}/igv.html'],
    source = [vardict_bgz,
              '${ref_dir}/${reference}.fa',
              mq_filtered_bam],
    action = ('$tabix -f ${SOURCES[0]}; '
	      '$igv ${SOURCES[0]} ${SOURCES[1]} --flanking $igv_flank '
              '--info-columns $igv_info '
              '--tracks ${SOURCES[0]} ${SOURCES[2]} --output ${TARGETS[1]}')
)

vardict_csv, vardict_bed = env.Command(
    target = ['$out/${sample}/${vardict_out}/sample_normalized.vcf.csv',
	      '$out/${sample}/${vardict_out}/sample_normalized.vcf.csv.bed'],
    source = vardict_normalized,
    action = ('$py_deps bin/to_csv.py $SOURCE ${TARGETS[0]} --sample $sample; '
	      '$py_deps bin/to_bed.py $TARGETS --split_mut $mutation_to_flank '
              '--bp $num_flanking_bp')
)

vardict_cov_bed, vardict_cov_csv = env.Command(
    target = ['$out/${sample}/${vardict_out}/'
              'sample_normalized_genomecov_intersect.bed',
              '$out/${sample}/${vardict_out}/sample_normalized_cov.csv'],
    source = [vardict_csv,
              vardict_bed,
              genome_cov],
    action = ('$bedtools intersect -a ${SOURCES[1]} -b ${SOURCES[2]} -wo '
              '> ${TARGETS[0]}; '
	      '$py_deps bin/add_cov.py ${TARGETS[0]} ${SOURCES[0]} '
              '${TARGETS[1]}')
)


# ################### Filter Variant Calls by DP ######################

gatk_cov_filtered = env.Command(
    target = '$out/${sample}/${gatk_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = gatk_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)

bcftools_cov_filtered = env.Command(
    target = '$out/${sample}/${bcftools_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = bcftools_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)

freebayes_cov_filtered = env.Command(
    target = '$out/${sample}/${freebayes_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = freebayes_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)

deepvariant_cov_filtered = env.Command(
    target = '$out/${sample}/${deepvariant_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = deepvariant_pass,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)

discosnp_cov_filtered = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = discosnp_pass_sorted,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)

lancet_cov_filtered = env.Command(
    target = '$out/${sample}/${lancet_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = lancet_final_vcf,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)

vardict_cov_filtered = env.Command(
    target = '$out/${sample}/${vardict_out}/'
             'sample_normalized_dp${min_read_depth}.vcf',
    source = vardict_normalized,
    action = '$bcftools bcftools filter -i \'FORMAT/DP>=${min_read_depth}\' '
             '-o $TARGET $SOURCE'
)


# ################### Check Calls ######################

gatk_calls = env.Command(
    target = '$out/${sample}/${gatk_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [gatk_cov_filtered,
              variant_cov_csv,
              gatk_cov_csv],
    action = '$py_deps bin/checker.py $gatk_out $sample $SOURCES $TARGET'
)

bcftools_calls = env.Command(
    target = '$out/${sample}/${bcftools_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [bcftools_cov_filtered,
              variant_cov_csv,
              bcftools_cov_csv],
    action = '$py_deps bin/checker.py $bcftools_out $sample $SOURCES $TARGET'
)

freebayes_calls = env.Command(
    target = '$out/${sample}/${freebayes_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [freebayes_cov_filtered,
              variant_cov_csv,
              freebayes_cov_csv],
    action = '$py_deps bin/checker.py $freebayes_out $sample $SOURCES $TARGET'
)

deepvariant_calls = env.Command(
    target = '$out/${sample}/${deepvariant_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [deepvariant_cov_filtered,
              variant_cov_csv,
              deepvariant_cov_csv],
    action = '$py_deps bin/checker.py $deepvariant_out $sample $SOURCES $TARGET'
)

discosnp_calls = env.Command(
    target = '$out/${sample}/${discosnp_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [discosnp_cov_filtered,
              variant_cov_csv,
              discosnp_cov_csv],
    action = '$py_deps bin/checker.py $discosnp_out $sample $SOURCES $TARGET'
)

lancet_calls = env.Command(
    target = '$out/${sample}/${lancet_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [lancet_cov_filtered,
              variant_cov_csv,
              lancet_cov_csv],
    action = '$py_deps bin/checker.py $lancet_out $sample $SOURCES $TARGET'
)

vardict_calls = env.Command(
    target = '$out/${sample}/${vardict_out}/'
             'sample_normalized_dp${min_read_depth}_checked.csv',
    source = [vardict_cov_filtered,
              variant_cov_csv,
              vardict_cov_csv],
    action = '$py_deps bin/checker.py $vardict_out $sample $SOURCES $TARGET'
)

all_checked_csv = env.Command(
    target = '$out/${sample}/sample_normalized_dp${min_read_depth}_checked.csv',
    source = [gatk_calls,
              bcftools_calls,
              freebayes_calls,
              deepvariant_calls,
              discosnp_calls,
              lancet_calls,
              vardict_calls],
    action = ('echo '
              '\'CHROM,POS,REF,ALT,TYPE,INS_TYPE,LEN,QUAL,AD_REF,AD_ALT,DP,'
              'BAM_DP,GT,ZYG,RK_DISCOSNP,TOOL,SAMPLE,TRUE_POS,FALSE_POS,'
              'FALSE_NEG\' > $TARGET; cat $SOURCES | sed '
              '\'/CHROM,POS,REF,ALT,TYPE,INS_TYPE,LEN,QUAL,AD_REF,AD_ALT,DP,'
              'BAM_DP,GT,ZYG,RK_DISCOSNP,TOOL,SAMPLE,TRUE_POS,FALSE_POS,'
              'FALSE_NEG/d\' >> $TARGET')
)

