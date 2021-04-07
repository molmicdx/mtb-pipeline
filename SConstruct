import os
import sys
from utils import parse_config

config = parse_config('configs/settings.conf')

# Ensure that a virtualenv is active before importing non-stdlib dependencies.
venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required'.format(venv))

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
vars.Add('log', '', 'logdir')
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
    reference=config.get('variant_simulation', 'reference_genome'),
    mutations_list=config.get('variant_simulation', 'mutations_list_output'),
    mutated_genome=config.get('variant_simulation', 'mutated_genome_output')
)

# Help(vars.GenerateHelpText(env))

# ############### start inputs ################
# ############# Simulate Variants #############
simulated_variants = env.Command(
    target=['$out/variants/GCF_000195955.2_20snps.txt', '$out/variants/GCF_000195955.2_20snps.fa'],
    source='$reference',
    action='python bin/variants.py $SOURCE $mutations_list $mutated_genome'
)

# ############### end inputs ##################
