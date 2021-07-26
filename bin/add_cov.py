import argparse
import csv

parser = argparse.ArgumentParser(description='Add read depth to true mutations file. Input bed file is produced by bedtools from the left outer join intersection of true variants bed and genome coverage bed files')

parser.add_argument('variant_bed', type=argparse.FileType('r'), help='variant bed file with genome coverage information')
parser.add_argument('true_mutations_csv', type=argparse.FileType('r'), help='csv file of true mutations')
parser.add_argument('variant_cov', type=argparse.FileType('w'), help='output csv with mutations and average read depths')
#parser.add_argument('cov_threshold', help='remove mutations with coverage below INT')

args = parser.parse_args()

true_mutations_reader = csv.DictReader(args.true_mutations_csv)
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

writer = csv.DictWriter(args.variant_cov, fieldnames=true_mutations_reader.fieldnames)
writer.writeheader()

true_mutation = next(true_mutations_reader, None)
while true_mutation:
    true_mutation['BAM_DP'] = round(sum(cov[int(true_mutation['POS'])])/len(cov[int(true_mutation['POS'])]))
    writer.writerow(true_mutation)
    true_mutation = next(true_mutations_reader, None)

