import argparse
import csv
import vcfpy
import sys

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (vcf), to csv for checker.py')
parser.add_argument('file', type=argparse.FileType('r'), help='vcf input file')
parser.add_argument('outcsv', type=argparse.FileType('w'), help='name of output csv file')

def vcf_to_csv(vcf, csvout):
    vcf_reader = vcfpy.Reader(vcf)
    fieldnames = ['CHROM','POS','REF','ALT','TYPE','INS_TYPE','AD_REF','AD_ALT','DP','BAM_DP','GT','ZYG','TRUE_POS','FALSE_POS','FALSE_NEG','TOOL','SAMPLE']
    writer = csv.DictWriter(csvout, fieldnames=fieldnames)
    writer.writeheader()
    record = next(vcf_reader, None)
    while record:
        true_variant = {}
        true_variant['CHROM'] = record.CHROM
        true_variant['POS'] = str(record.POS)
        true_variant['REF'] = record.REF
        true_variant['ALT'] = ','.join([alt.value for alt in record.ALT])
        true_variant['TYPE'] = record.INFO['TYPE'][0]
        try:
            true_variant['INS_TYPE'] = record.INFO['INS_TYPE'][0]
        except KeyError:
            true_variant['INS_TYPE'] = None
        for call in record.calls:
            true_variant['SAMPLE'] = call.sample
            true_variant['GT'] = call.data['GT']
            true_variant['ZYG'] = 'hom'
            true_variant['TRUE_POS'] = 0
            true_variant['FALSE_POS'] = 0
            true_variant['FALSE_NEG'] = 1
        writer.writerow(true_variant)
        record = next(vcf_reader, None)
    return

def main():
    args = parser.parse_args()
    vcf_to_csv(args.file, args.outcsv)

if __name__ == '__main__':
    sys.exit(main())
