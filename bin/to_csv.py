import argparse
import csv
import vcfpy

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (vcf), to csv for checker.py')

parser.add_argument('file', type=argparse.FileType('r'), help='vcf input file')
parser.add_argument('outcsv', type=argparse.FileType('w'), help='name of output csv file')
parser.add_argument('-t', '--txt', action='store_true', help='input is text file')
args = parser.parse_args()

vcf_reader = vcfpy.Reader(args.file)
fieldnames = ['CHROM','POS','REF','ALT','TYPE','AD_REF','AD_ALT','DP','BAM_DP','GT','ZYG','TRUE_POS','FALSE_POS','FALSE_NEG','TOOL','SAMPLE']
writer = csv.DictWriter(args.outcsv, fieldnames=fieldnames)
writer.writeheader()
record = next(vcf_reader, None)
while record:
    true_variant = {}
    true_variant['CHROM'] = record.CHROM
    true_variant['POS'] = str(record.POS)
    true_variant['REF'] = record.REF
    true_variant['ALT'] = ','.join([alt.value for alt in record.ALT])
    for call in record.calls:
        true_variant['SAMPLE'] = call.sample
        true_variant['GT'] = call.data['GT']
        true_variant['ZYG'] = 'hom'
        true_variant['TRUE_POS'] = 0
        true_variant['FALSE_POS'] = 0
        true_variant['FALSE_NEG'] = 1
    if len(true_variant['REF']) == 1 and len(true_variant['ALT']) == 1:
        true_variant['TYPE'] = 'SNP'
    elif len(true_variant['REF']) > len(true_variant['ALT']):
        true_variant['TYPE'] = 'DEL'
    elif len(true_variant['REF']) < len(true_variant['ALT']):
        true_variant['TYPE'] = 'INS'
    writer.writerow(true_variant)
    record = next(vcf_reader, None)
