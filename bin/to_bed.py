import argparse

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (csv) to bed file')

parser.add_argument('file', help='csv of variants.py mutation list')
args = parser.parse_args()

bedfile = args.file + '.bed'
with open(args.file, 'r') as mutation_csv:
    mutations = mutation_csv.readlines()
    with open(bedfile, 'w') as outbed:
        for row in mutations[1:]:
            data = row.split(',')
            chromStart = int(data[1]) - 1
            size = abs(len(data[3].rstrip()) - len(data[2].rstrip())) #length of mutation
            chromEnd = chromStart + size
            outbed.write(data[0] + '\t' + str(chromStart) + '\t' + str(chromEnd) + '\n')

