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
            pos = int(data[1])
            size = len(data[3].rstrip())
            chromEnd = str(pos + size)
            outbed.write('\t'.join(data[:2]) + '\t' + chromEnd + '\n')

