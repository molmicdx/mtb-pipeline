import os
import sys
from utils import parse_config

config = parse_config('configs/settings.conf')

# Ensure that a virtualenv is active before importing non-stdlib dependencies.
venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required'.format(venv))

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
    ENV=dict(os.environ, PATH=PATH, SHELLOPTS='errexit:pipefail'),
    variables=vars,
    SHELL='bash',
    cwd=os.getcwd(),
    reference=config.get('DEFAULT', 'reference_genome'),
    variant=config.get('variant_simulation', 'variant_name'),
    variants_config=config.get('variant_simulation', 'variants_config'),
    variants_out=config.get('output', 'variants'),
    gatk='{} exec -B $cwd {} gatk'.format(singularity, gatk_img)
)

# Help(vars.GenerateHelpText(env))

# ############### start inputs ################
# ############# Simulate Variants #############
simulated_variants_table, simulated_variants_fa = env.Command(
    target=['$out/$variants_out/${variant}.txt', '$out/$variants_out/${variant}.fa'],
    source='$reference',
    action=('python bin/variants.py --settings $variants_config $SOURCE $TARGETS')
)

# ############# Normalize Variant VCF ###############
simulated_variant_vcf = env.Command(
    target='$out/$variants_out/${variant}.txt.vcf',
    source=simulated_variants_table,
    action='python bin/to_vcf.py $SOURCE $variant'
)

normalized_variant_vcf = env.Command(
    target=['$out/$variants_out/${variant}_normalized.vcf', '$log/$variants_out/${variant}_normalized.log'],
    source=simulated_variant_vcf,
    action=('$gatk LeftAlignAndTrimVariants -R $reference -V $SOURCE -O ${TARGETS[0]} '
            '> ${TARGETS[1]} 2<&1')
)

# ############### end inputs ##################
