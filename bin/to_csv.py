import argparse

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (txt) or any vcf, to csv for checker.py')

parser.add_argument('file', help='txt or vcf input file')
parser.add_argument('sample', help='sample name')
parser.add_argument('-t', '--txt', action='store_true', help='input is txt')
parser.add_argument('-v', '--vcf', action='store_true', help='input is vcf')
args = parser.parse_args()

csvfile = args.file + '.csv'
csv_header = 'CHROM,POS,REF,ALT,' + args.sample + '\n'

with open(args.file, 'r') as mutation_list:
    mutations = mutation_list.readlines()
    with open(csvfile, 'w') as outcsv:
        outcsv.write(csv_header)
        if args.txt:
            for row in mutations[1:]:
                # check that it is not an empty row
                if len(row) > 1:
                    csv_row_fields = row.rstrip().split('\t')
                    csv_row_fields.pop(2)
                    csv_row = ','.join(csv_row_fields) + ',1\n'
                    outcsv.write(csv_row)
        elif args.vcf:
            vcf_mutations = []
            for line in mutations:
                if line.startswith('#'):
                    continue
                else:
                    vcf_mutations.append(line)
            for row in vcf_mutations:
                # check that it is not an empty row
                if len(row) > 1:
                    csv_row_fields = row.rstrip().split('\t')
                    new_csv_row = csv_row_fields[:2] + csv_row_fields[3:5]
                    csv_row = ','.join(new_csv_row) + ',1\n'
                    outcsv.write(csv_row)

