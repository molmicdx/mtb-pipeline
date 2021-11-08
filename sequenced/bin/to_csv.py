import argparse
import csv
import vcfpy
import sys

parser = argparse.ArgumentParser(description='Convert variants.py mutation list (vcf) or csv mutation list (CHROM,POS,REF,ALT,TYPE,<sample>), to csv for checker.py')
parser.add_argument('file', type=argparse.FileType('r'), help='input file')
parser.add_argument('outcsv', type=argparse.FileType('w'), help='name of output csv file')
parser.add_argument('--formatcsv', action='store_const', const=True, help='input file is csv')
parser.add_argument('--asref', action='store_const', const=True,  help='use with --formatcsv if changing introduced mutations to true mutations csv')
parser.add_argument('--sample', default=None, help='provide new sample name if using --asref')

def format_csv(csvin, csvout, toref, samplename):
    reader = csv.DictReader(csvin)
    fields_in = reader.fieldnames
    fieldnames = ['CHROM','POS','REF','ALT','TYPE','INS_TYPE','LEN','QUAL','AD_REF','AD_ALT','DP','BAM_DP','GT','ZYG','RK_DISCOSNP','TOOL','SAMPLE','TRUE_POS','FALSE_POS','FALSE_NEG']
    writer = csv.DictWriter(csvout, fieldnames=fieldnames)
    writer.writeheader()
    mutation = next(reader, None)
    while mutation:
        entry = {}
        entry['CHROM'] = mutation['CHROM']
        entry['POS'] = mutation['POS']
        if toref:
            entry['REF'] = mutation['ALT']
            entry['ALT'] = mutation['REF']
            if mutation['TYPE'] != 'SNP':
                if mutation['TYPE'] == 'DEL':
                    entry['TYPE'] = 'INS'
                elif mutation['TYPE'] == 'INS':
                    entry['TYPE'] = 'DEL'
                    entry['INS_TYPE'] = ''
            else:
                entry['TYPE'] = mutation['TYPE']
        else:
            entry['REF'] = mutation['REF']
            entry['ALT'] = mutation['ALT']
            entry['TYPE'] = mutation['TYPE']
        entry['GT'] = mutation[fields_in[-1]]
        entry['ZYG'] = 'hom'
        entry['TRUE_POS'] = 0
        entry['FALSE_POS'] = 0
        entry['FALSE_NEG'] = 1
        if samplename == None:
            try:
                entry['SAMPLE'] = mutation['SAMPLE']
            except KeyError:
                entry['SAMPLE'] = fields_in[-1]
        else:
            entry['SAMPLE'] = samplename
        writer.writerow(entry)
        mutation = next(reader, None)
    
    return

def vcf_to_csv(vcf, csvout):
    vcf_reader = vcfpy.Reader(vcf)
    fieldnames = ['CHROM','POS','REF','ALT','TYPE','INS_TYPE','LEN','QUAL','AD_REF','AD_ALT','DP','BAM_DP','GT','ZYG','RK_DISCOSNP','TOOL','SAMPLE','TRUE_POS','FALSE_POS','FALSE_NEG']
    writer = csv.DictWriter(csvout, fieldnames=fieldnames)
    writer.writeheader()
    record = next(vcf_reader, None)
    while record:
        variant = {}
        variant['CHROM'] = record.CHROM
        variant['POS'] = str(record.POS)
        variant['REF'] = record.REF
        variant['ALT'] = ','.join([alt.value for alt in record.ALT])
        try:
            variant['TYPE'] = ''.join(record.INFO['TYPE'])[:3].upper() # attempt to standardize TYPEs reported by different variant callers
        except KeyError:
            size = len(variant['REF'].rstrip()) - len(variant['ALT'].rstrip())
            if size == 0: # SNP size
                variant['TYPE'] = 'SNP'
            elif size > 0:
                variant['TYPE'] = 'DEL'
            elif size < 0:
                variant['TYPE'] = 'INS'
        try:
            variant['INS_TYPE'] = record.INFO['INS_TYPE'][0]
        except KeyError:
            variant['INS_TYPE'] = None
        for call in record.calls:
            variant['SAMPLE'] = call.sample
            variant['GT'] = str(call.data['GT'])
            if len(variant['GT']) > 1:
                if variant['GT'][0] != variant['GT'][-1]:
                    variant['ZYG'] = 'het'
                else:
                    variant['ZYG'] = 'hom'
            variant['TRUE_POS'] = 0
            variant['FALSE_POS'] = 0
            variant['FALSE_NEG'] = 1
        writer.writerow(variant)
        record = next(vcf_reader, None)
    return

def main():
    args = parser.parse_args()
    if args.formatcsv:
        format_csv(args.file, args.outcsv, args.asref, args.sample)
    else:
        vcf_to_csv(args.file, args.outcsv)

if __name__ == '__main__':
    sys.exit(main())
