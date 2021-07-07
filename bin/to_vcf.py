import argparse

parser = argparse.ArgumentParser(description='Convert variants.py mutation list to vcf tolerated by gatk v4.0.11.0')

parser.add_argument('file', help='variants.py txt file output of introduced mutations')
parser.add_argument('sample', help='sample name')
parser.add_argument('-c', '--csv', action='store_true', help='output csv for checker.py')
args = parser.parse_args()

with open(args.file, 'r') as mutation_list:
    header_fields = mutation_list.readline().split('\t')
    header = '\t'.join(header_fields[:-1])
    new_header = header.upper().rstrip() + '\tQUAL\tFILTER\tINFO\tFORMAT\t' + args.sample + '\n'
    if new_header[0] != '#':
        new_header = '#' + new_header        
    mutations = mutation_list.readlines()

outfile = args.file + '.vcf'
with open(outfile, 'w') as outvcf:
    outvcf.write('##fileformat=VCFv4.1\n##INFO=<ID=TY,Number=A,Type=String,Description="SNP, INS, or DEL">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    outvcf.write(new_header)
    for row in mutations:
        # check that it is not an empty row
        if len(row) > 1:
            first_5_cols = row.rstrip().split('\t')[:5]
            TYPE = row.rstrip().split('\t')[-1]
            new_row = '\t'.join(first_5_cols) + '\t.\t.\tTY=' + TYPE + '\tGT\t1\n'
            outvcf.write(new_row)

if args.csv:
    csvfile = args.file + '.csv'
    csv_header = 'CHROM,POS,REF,ALT,TYPE,' + args.sample + '\n'
    with open(csvfile, 'w') as outcsv:
        outcsv.write(csv_header)
        for row in mutations:
            # check that it is not an empty row
            if len(row) > 1:
                csv_row_fields = row.rstrip().split('\t')
                csv_row_fields.pop(2)
                csv_row = ','.join(csv_row_fields) + ',1\n'
            outcsv.write(csv_row)

