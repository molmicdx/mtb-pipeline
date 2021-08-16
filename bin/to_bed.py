import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (csv) to bed file')
parser.add_argument('file', type=argparse.FileType('r'), help='csv of variants.py mutation list')
parser.add_argument('outbed', type=argparse.FileType('w'), help='name of output bed file')
parser.add_argument('--split_mut', default='DEL', help='type of mutation to get flanking read depths of')
parser.add_argument('--bp', default=5, help='number of flanking base pairs to get read depths for type of mutation specified in split_mut')

def get_variant_bed(mutation):
    variant_bed = {}
    variant_bed['chrom'] = mutation['CHROM']
    variant_bed['type'] = mutation['TYPE']
    variant_bed['chromStart'] = int(mutation['POS']) - 1
    size = abs(len(mutation['ALT'].rstrip()) - len(mutation['REF'].rstrip()))
    if size == 0: # SNP size
        size = 1
    variant_bed['chromEnd'] = variant_bed['chromStart'] + size + 1
    return variant_bed

def split_call(variant_bed, bp):
    variant_bed1 = variant_bed.copy()
    variant_bed1['chromEnd'] = variant_bed['chromStart']
    variant_bed1['chromStart'] -= int(bp)
    variant_bed2 = variant_bed.copy()
    variant_bed2['chromStart'] = variant_bed['chromEnd']
    variant_bed2['chromEnd'] += int(bp)

    return variant_bed1, variant_bed2


def write_to_bed(csv_in, bed_out, muttype, bp):
    all_mutations = csv.DictReader(csv_in)
    fieldnames = ['chrom', 'chromStart', 'chromEnd', 'genomeCov', 'type']
    bed_writer = csv.DictWriter(bed_out, fieldnames=fieldnames, delimiter='\t')
    mutation = next(all_mutations, None)
    while mutation:
        variant_call = get_variant_bed(mutation)
        if variant_call['type'] == muttype:
            before_mut, after_mut = split_call(variant_call, bp)
            bed_writer.writerow(before_mut)
            bed_writer.writerow(variant_call)
            bed_writer.writerow(after_mut)
        else:
            bed_writer.writerow(variant_call)
        mutation = next(all_mutations, None)            
    return

def main():
    args = parser.parse_args()
    write_to_bed(args.file, args.outbed, args.split_mut, args.bp)

if __name__ == '__main__':
    sys.exit(main())

