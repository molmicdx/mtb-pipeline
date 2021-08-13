import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (csv) to bed file')
parser.add_argument('file', type=argparse.FileType('r'), help='csv of variants.py mutation list')
parser.add_argument('outbed', type=argparse.FileType('w'), help='name of output bed file')

def get_variant_bed(mutation):
    variant_bed = {}
    variant_bed['chrom'] = mutation['CHROM']
    variant_bed['chromStart'] = int(mutation['POS']) - 1
    size = abs(len(mutation['ALT'].rstrip()) - len(mutation['REF'].rstrip()))
    if size == 0: # SNP size
        size = 1
    variant_bed['chromEnd'] = variant_bed['chromStart'] + size + 1
    return variant_bed


def write_to_bed(csv_in, bed_out):
    all_mutations = csv.DictReader(csv_in)
    fieldnames = ['chrom', 'chromStart', 'chromEnd']
    bed_writer = csv.DictWriter(bed_out, fieldnames=fieldnames, delimiter='\t')
    mutation = next(all_mutations, None)
    while mutation:
        bed_writer.writerow(get_variant_bed(mutation))
        mutation = next(all_mutations, None)            
    return

def main():
    args = parser.parse_args()
    write_to_bed(args.file, args.outbed)

if __name__ == '__main__':
    sys.exit(main())


'''
args = parser.parse_args()

bedfile = args.file + '.bed'
with open(args.file, 'r') as mutation_csv:
    mutations = mutation_csv.readlines()
    with open(bedfile, 'w') as outbed:
        for row in mutations[1:]:
            data = row.split(',')
            chromStart = int(data[1]) - 1 #BED coordinates are 0-start, half-open
            size = abs(len(data[3].rstrip()) - len(data[2].rstrip())) #length of mutation
            if size == 0: #SNP size
                size = 1
            chromEnd = chromStart + size + 1
            outbed.write(data[0] + '\t' + str(chromStart) + '\t' + str(chromEnd) + '\n')
'''
