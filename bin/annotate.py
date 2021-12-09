from binary_tree import *
import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='annotate with genome loci')
parser.add_argument('input_csv', type=argparse.FileType('r'), help='final checked csv')
parser.add_argument('loci_list_csv', help='csv list of locus,start,end')
parser.add_argument('output_csv', type=argparse.FileType('w'), help='output csv with locus annotation')

def main():
    args = parser.parse_args()
    with open(args.loci_list_csv, 'r') as amr:
        lines = amr.readlines()

    features = []
    for each in lines:
        feat = each.rstrip().split(',')
        # Convert BED coordinates to VCF coordinates
        features.append((feat[0], int(feat[1])+1, int(feat[2]))) 
    mytree = RangeTree(features)

    variants_reader = csv.DictReader(args.input_csv)
    fieldnames = variants_reader.fieldnames
    fieldnames.append('LOCI')
    writer = csv.DictWriter(args.output_csv, fieldnames = fieldnames)
    writer.writeheader()
    variant = next(variants_reader, None)
    while variant:
        variant['LOCI'] = mytree.assign(int(variant['POS']))
        writer.writerow(variant)
        variant = next(variants_reader, None)
    #print(mytree.assign(5245))

if __name__ == '__main__':
    sys.exit(main())


