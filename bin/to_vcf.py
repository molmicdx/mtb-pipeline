import argparse

parser = argparse.ArgumentParser(description='Convert variants.py mutation list to vcf tolerated by gatk v4.0.11.0')

parser.add_argument('file', help='variants.py txt file output of introduced mutations')
parser.add_argument('sample', help='sample name')
args = parser.parse_args()

with open(args.file, 'r') as mutation_list:
    header = mutation_list.readline()
    new_header = header.upper().rstrip() + '\tQUAL\tFILTER\tINFO\tFORMAT\t' + args.sample + '\n'
    if new_header[0] != '#':
        new_header = '#' + new_header
    mutations = mutation_list.readlines()

outfile = args.file + '.vcf'
with open(outfile, 'w') as outvcf:
    outvcf.write('##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    outvcf.write(new_header)
    for row in mutations:
        # check that it is not an empty row
        if len(row) > 1:
            new_row = row.rstrip() + '\t.\t.\t.\tGT\t1\n'
            outvcf.write(new_row)

