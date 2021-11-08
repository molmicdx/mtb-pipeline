from configparser import SafeConfigParser
import csv
import os
import sys


def parse_config(filename):
    """ Parse config file, exit if it can't be found """
    if not os.path.exists(filename):
        sys.exit("Cannot find configuration file '{}'".format(filename))
    config = SafeConfigParser(allow_no_value=True)
    config.optionxform = str # Make sure keys are case sensitive
    config.read(filename)
    config.read(config.get('sequenced', 'reads_data'))
    return config


def annotate_specimens(config):
    """ Retrieve list of specimens and their metadata """
    specimens = {}
    for row in csv.DictReader(open(config.get('sequenced', 'manifest'), 'r')):
        specimens["{}-{}".format(row['run_number'], row['barcode_id'])] = row
    for specimen_id, metadata in specimens.items():
        metadata['specimen_id'] = specimen_id
        if not metadata['label']:
            metadata['label'] = specimen_id
        reads = dict(config.items(specimen_id))
        metadata['forward_reads'] = reads['forward']
        metadata['reverse_reads'] = reads['reverse']
        if not os.path.exists(reads['forward']) or not os.path.exists(reads['reverse']):
            sys.exit("Files with reads for specimen {} don't exist".format(specimen_id))
    return specimens


def write_matrix_to_csv(names, matrix, filename):
    """ Writes matrix as csv with names as top row and leftmost column """
    with open(filename, 'w') as out:
        out.write('Sample,' + ','.join(names) + '\n')
        for name1 in names:
            out.write(name1 + ',' + ','.join([str(matrix[name1][name2]) for name2 in names]) + '\n')

#def main():
#    config = parse_config('/molmicro/working/ymseah/mtb_amr/sequenced/configs/settings.conf')
#    print(annotate_specimens(config))

#if __name__=='__main__':
#    sys.exit(main())

