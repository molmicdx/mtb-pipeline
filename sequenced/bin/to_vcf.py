import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Convert variants.py mutation list to vcf tolerated by gatk v4.0.11.0')
parser.add_argument('file', help='variants.py txt file output of introduced mutations')
parser.add_argument('sample', help='sample name')


def get_variant(mutation):
    args = parser.parse_args()
    variant = {}
    variant['#CHROM'] = mutation['#CHROM']
    variant['POS'] = mutation['POS']
    variant['ID'] = '.'
    variant['QUAL'] = '.'
    variant['FILTER'] = '.'
    variant['REF'] = mutation['REF']
    variant['ALT'] = mutation['ALT']
    variant['FORMAT'] = 'GT'
    variant[args.sample] = '1'
    try:
        variant['INFO'] = 'TYPE=' + mutation['TYPE'] + ';INS_TYPE=' + mutation['INS_TYPE']
    except TypeError:
        variant['INFO'] = 'TYPE=' + mutation['TYPE']
    return variant

def main():
    args = parser.parse_args()
    with open(args.file, 'r', newline='') as infile:
        mutations_reader = csv.DictReader(infile, delimiter='\t')
        vcf_header = '##fileformat=VCFv4.1\n##INFO=<ID=TYPE,Number=A,Type=String,Description="SNP, INS, or DEL">\n##INFO=<ID=INS_TYPE,Number=A,Type=String,Description="DUP, INV, or RDM">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##source=to_vcf.py\n'
        with open(args.file + '.vcf', 'w', newline='') as outvcf:
            outvcf.write(vcf_header)
            fieldnames = mutations_reader.fieldnames[:-2] + ['QUAL', 'FILTER', 'INFO', 'FORMAT', args.sample]
            variant_writer = csv.DictWriter(outvcf, fieldnames=fieldnames, delimiter='\t')
            variant_writer.writeheader()
            mutation = next(mutations_reader, None)
            while mutation:
                variant_writer.writerow(get_variant(mutation))
                mutation = next(mutations_reader, None)

if __name__ == '__main__':
    sys.exit(main())

