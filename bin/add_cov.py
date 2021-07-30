import argparse
import csv

parser = argparse.ArgumentParser(description='Add read depth to variants csv file. Input bed file is produced by bedtools from the left outer join intersection of variants bed and genome coverage bed files')

parser.add_argument('variant_bed', type=argparse.FileType('r'), help='variant bed file with genome coverage information')
parser.add_argument('variants_csv', type=argparse.FileType('r'), help='csv file of variants')
parser.add_argument('variant_cov', type=argparse.FileType('w'), help='output csv with variants and average read depths (i.e. BAM_DP)')
parser.add_argument('cov_threshold', help='remove variants with coverage below INT')

args = parser.parse_args()

variants_reader = csv.DictReader(args.variants_csv)
mutation_cov = next(csv.reader(args.variant_bed, delimiter='\t'), None)
cov = {}
while mutation_cov:
    # cov dictionary key is (chromStart + 1): value is list of read depths.
    # BED uses zero-based positions, VCF uses 1-based positions
    try:
        cov[int(mutation_cov[1]) + 1].append(int(mutation_cov[-1]))
    except KeyError:
        cov[int(mutation_cov[1]) + 1] = [int(mutation_cov[-1])]
    mutation_cov = next(csv.reader(args.variant_bed, delimiter='\t'), None)

writer = csv.DictWriter(args.variant_cov, fieldnames=variants_reader.fieldnames)
writer.writeheader()

variant = next(variants_reader, None)
while variant:
    variant['BAM_DP'] = round(sum(cov[int(variant['POS'])])/len(cov[int(variant['POS'])]))
    if variant['BAM_DP'] >= int(args.cov_threshold):
        writer.writerow(variant)
    variant = next(variants_reader, None)

