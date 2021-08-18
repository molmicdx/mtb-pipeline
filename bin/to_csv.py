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
        variant = {}
        variant['CHROM'] = record.CHROM
        variant['POS'] = str(record.POS)
        variant['REF'] = record.REF
        variant['ALT'] = ','.join([alt.value for alt in record.ALT])
        try:
            variant['TYPE'] = ''.join(record.INFO['TYPE'])
        except KeyError:
            size = len(variant['REF'].rstrip()) - len(variant['ALT'].rstrip())
            if size == 0: # SNP size
                variant['TYPE'] = 'SNP'
            elif size > 0:
                variant['TYPE'] = 'DEL'
            elif size < 0:
                variant['TYPE'] = 'INS'
        try:
            variant['INS_TYPE'] = record.INFO['INS_TYPE'][0]
        except KeyError:
            variant['INS_TYPE'] = None
        for call in record.calls:
            variant['SAMPLE'] = call.sample
            variant['GT'] = call.data['GT']
            variant['ZYG'] = 'hom'
            variant['TRUE_POS'] = 0
            variant['FALSE_POS'] = 0
            variant['FALSE_NEG'] = 1
        writer.writerow(variant)
        record = next(vcf_reader, None)
    return

def main():
    args = parser.parse_args()
    vcf_to_csv(args.file, args.outcsv)

if __name__ == '__main__':
    sys.exit(main())
