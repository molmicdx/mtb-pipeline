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
    variant = config.get('variant_simulation', 'variant_name'),
    variants_config = config.get('variant_simulation', 'variants_config'),
    read1_q = config.get('variant_simulation', 'r1_qual_profile'),
    read2_q = config.get('variant_simulation', 'r2_qual_profile'),
    read_len = config.get('variant_simulation', 'read_length'),
    read_depth = config.get('variant_simulation', 'avg_read_depth'),
    frag_len = config.get('variant_simulation', 'avg_fragment_length'),
    sd = config.get('variant_simulation', 'std_dev'),
    variants_out = config.get('output', 'variants'),
    reads_out = config.get('output', 'reads'),
    gatk = '{} exec -B $cwd {} gatk'.format(singularity, gatk_img)
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
    action = ('gzip $SOURCES')
)




# ############### end inputs ##################
