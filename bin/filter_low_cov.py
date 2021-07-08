import argparse

parser = argparse.ArgumentParser(description='Filter true mutations based on read depth')

parser.add_argument('variant_bed', help='variant bed file with coverage information')
parser.add_argument('true_mutations_csv', help='csv file of true mutations')
parser.add_argument('cov_threshold', help='remove mutations with coverage below INT')
args = parser.parse_args()

cov_dict={}
with open(args.variant_bed, 'r') as covfile:
    bedlines = covfile.readlines()
    for line in bedlines:
        # Create dictionary of position:coverage
        cov_dict[line.split('\t')[1]] = int(line.split('\t')[-1])

filtered_file = args.true_mutations_csv.split('.')[0] + '_dp'+ args.cov_threshold + '.csv'
with open(args.true_mutations_csv, 'r') as mutation_list:
    mutations = mutation_list.readlines()
    with open(filtered_file, 'w') as outcsv:
        outcsv.write(mutations[0])
        for row in mutations[1:]:
            if row.split(',')[4] == 'DEL':
                outcsv.write(row)
            else:
                if cov_dict[row.split(',')[1]] > int(args.cov_threshold):
                    outcsv.write(row)

