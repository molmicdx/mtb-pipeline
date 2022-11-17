#!/usr/bin/env python3

from configparser import ConfigParser
from datetime import datetime
import os
import csv
import sys
import re
import subprocess

def parse_config(filename):
    """ Parse config file, exit if it can't be found """
    if not os.path.exists(filename):
        sys.exit("Cannot find configuration file '{}'".format(filename))
    config = ConfigParser(allow_no_value=True)
    config.optionxform = str # Make sure keys are case sensitive
    config.read(filename)
    #try: 
    #    config.read(config.get('sequenced', 'reads_data'))
    #except Exception:
    #    print("No sequenced reads. Sequences will be simulated.")
    return config


def annotate_specimens(config):
    """ Retrieve list of specimens and their metadata """
    with open(config.get('sequenced', 'reads_data'), 'r') as f:
        fastq_list = f.read().splitlines()

    specimens = {}
    for row in csv.DictReader(open(config.get('sequenced', 'manifest'), 'r')):
        specimens["{}-{}".format(row['run_number'], row['barcode_id'])] = row
    for specimen_id, metadata in specimens.items():
        # define regex patterns for matching each of fwd and rev reads
        fwd_regex = re.compile(f".*{metadata['run_number']}[_-]{metadata['barcode_id']}_.*_R1.*.fastq.gz")
        rev_regex = re.compile(f".*{metadata['run_number']}[_-]{metadata['barcode_id']}_.*_R2.*.fastq.gz") 
        metadata['specimen_id'] = specimen_id
        if not metadata['label']:
            metadata['label'] = specimen_id
        fwd_reads = list(filter(fwd_regex.match, fastq_list))
        rev_reads = list(filter(rev_regex.match, fastq_list))
        # assert only one read file for each direction
        assert len(fwd_reads) == 1 and len(rev_reads) == 1

        metadata['forward'] = fwd_reads[0]
        metadata['reverse'] = rev_reads[0]
        if not os.path.exists(metadata['forward']) or not os.path.exists(metadata['reverse']):
            sys.exit("Files with reads for specimen {} don't exist".format(specimen_id))
    return specimens

def timestamp_now():
    return datetime.now().strftime("%A, %B, %d, %Y, %I:%M %p")

def get_versions(config_name):
    config = parse_config(config_name)
    versions = {}
    versions['pipeline'] = subprocess.check_output(\
            ["git", "describe", "--tags", "--dirty"]\
            ).strip().decode('UTF-8')
    versions['who'] = config.get('DEFAULT', 'reference_db')
    return versions
