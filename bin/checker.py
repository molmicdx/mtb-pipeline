import argparse
import csv
import sys
import vcfpy

def get_args():
    parser = argparse.ArgumentParser(description="Check pipeline results with true variants.")
    parser.add_argument('variant_caller', help="name of variant caller")
    parser.add_argument('sample', help="name of sample")
    parser.add_argument('merged_vcf', type=argparse.FileType('r'),
                        help="merged vcf file with pipeline results")
    parser.add_argument('true_variants', type=argparse.FileType('r'),
                        help="csv file with sorted true variants")
    parser.add_argument('variant_cov', type=argparse.FileType('r'),
                        help='csv file of called variants with BAM read depths')
    parser.add_argument('called_variants', type=argparse.FileType('w'),
                        help="output csv with called variants")
    #parser.add_argument('summary', help="output summary csv")
    return parser.parse_args()


def get_variant_from_vcf_record(record):
    args = get_args()
    variant = {}
    variant['CHROM'] = record.CHROM
    variant['POS'] = str(record.POS)
    variant['REF'] = record.REF
    variant['ALT'] = ','.join([alt.value for alt in record.ALT])
    variant['QUAL'] = record.QUAL
    try:
        variant['RK_DISCOSNP'] = record.INFO['Rk']
    except KeyError:
        variant['RK_DISCOSNP'] = ''
    for call in record.calls:
        variant['SAMPLE'] = call.sample
        variant['GT'] = str(call.data['GT'])
        if len(variant['GT']) > 1:
            if variant['GT'][0] != variant['GT'][-1]:
                variant['ZYG'] = 'het'
            else:
                variant['ZYG'] = 'hom'
        elif len(variant['GT']) == 1:
            variant['ZYG'] = 'hom'
        try:
            variant['DP'] = str(call.data['DP'])
        except KeyError: # bcftools doesn't have DP yet
            variant['DP'] = ''
        try:
            variant['AD_REF'] = str(call.data['AD'][0])
            variant['AD_ALT'] = str(call.data['AD'][1])
        except KeyError: #bcftools doesn't have AD yet
            variant['AD_REF'] = ''
            variant['AD_ALT'] = ''
        variant['TRUE_POS'] = 0
        variant['FALSE_POS'] = 0
        variant['FALSE_NEG'] = 0
        variant['TOOL'] = args.variant_caller
    
    return variant


def is_match(vcf_record, true_variant):
    chrom = vcf_record.CHROM == true_variant['CHROM']
    pos = vcf_record.POS == int(true_variant['POS'])
    ref = vcf_record.REF == true_variant['REF']
    alt = vcf_record.ALT[0].value == true_variant['ALT']
    #type = ','.join(vcf_record.INFO['TY']) == true_variant['TYPE']
    genotypes = True
    for call in vcf_record.calls:
        try: #Catch diploid GT reporting
            if call.data['GT'][0] != true_variant['GT'] and call.data['GT'][2] != true_variant['GT']:
                genotypes = False
        except IndexError:
            if call.data['GT'][0] and call.data['GT'][0] != true_variant['GT']:
                genotypes = False
    return chrom and pos and ref and alt and genotypes


def check(vcf_reader, true_variants_reader, vcf_cov_reader):
    args = get_args()
    all_variants = []
    tps, fps, fns = [], [], []
    vcf_record = next(vcf_reader, None)
    true_variant = next(true_variants_reader, None)
    vcf_cov = next(vcf_cov_reader, None)
    while vcf_record and true_variant:
        called = get_variant_from_vcf_record(vcf_record)
        size = len(called['REF']) - len(called['ALT'])
        if size == 0: 
            called['LEN'] = 1 # SNP size
        else:
            called['LEN'] = abs(size)
        while called['POS'] != vcf_cov['POS']:
            vcf_cov = next(vcf_cov_reader, None)
        if called['POS'] == vcf_cov['POS']:
            called['BAM_DP'] = vcf_cov['BAM_DP']
        if (vcf_record.CHROM, vcf_record.POS) <= (true_variant['CHROM'], int(true_variant['POS'])):
            if is_match(vcf_record, true_variant):
                called['TRUE_POS'] = 1
                called['TYPE'] = true_variant['TYPE']
                if called['TYPE'] == 'INS':
                    called['INS_TYPE'] = true_variant['INS_TYPE']
                all_variants.append(called)
                tps.append(called)
                true_variant = next(true_variants_reader, None)
            else:
                called['FALSE_POS'] = 1 # need fix for case when variant present in VCF but GT=0; filter out GT=0?
                if size == 0:
                    called['TYPE'] = 'SNP'
                    called['LEN'] = 1
                elif size > 0:
                    called['TYPE'] = 'DEL'
                    called['LEN'] = size
                elif size < 0:
                    called['TYPE'] = 'INS'
                    called['LEN'] = abs(size)
                all_variants.append(called)
                fps.append(called)
            vcf_record = next(vcf_reader, None)
            vcf_cov = next(vcf_cov_reader, None)
        else:
            true_variant['TOOL'] = args.variant_caller
            if abs(len(true_variant['REF']) - len(true_variant['ALT'])) == 0:
                true_variant['LEN'] = 1 # SNP size
            else:
                true_variant['LEN'] = abs(len(true_variant['REF']) - len(true_variant['ALT']))
            all_variants.append(true_variant)
            fns.append(true_variant)
            true_variant = next(true_variants_reader, None)
    while vcf_record:
        called = get_variant_from_vcf_record(vcf_record)
        size = len(called['REF']) - len(called['ALT'])
        while called['POS'] != vcf_cov['POS']:
            vcf_cov = next(vcf_cov_reader, None)
        if called['POS'] == vcf_cov['POS']:
            called['BAM_DP'] = vcf_cov['BAM_DP']
        called['FALSE_POS'] = 1
        if size == 0:
            called['TYPE'] = 'SNP'
            called['LEN'] = 1
        elif size > 0:
            called['TYPE'] = 'DEL'
            called['LEN'] = size
        elif size < 0:
            called['TYPE'] = 'INS'
            called['LEN'] = abs(size)
        all_variants.append(called)
        fps.append(called)
        vcf_record = next(vcf_reader, None)
        vcf_cov = next(vcf_cov_reader, None)
    while true_variant:
        true_variant['TOOL'] = args.variant_caller
        if abs(len(true_variant['REF']) - len(true_variant['ALT'])) == 0:
            true_variant['LEN'] = 1
        else:
            true_variant['LEN'] = abs(len(true_variant['REF']) - len(true_variant['ALT']))
        all_variants.append(true_variant)
        fns.append(true_variant)
        true_variant = next(true_variants_reader, None)
    
    return all_variants, tps, fps, fns


def write_variants(variants, fieldnames, file):
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(variants)


def main():
    args = get_args()
    all_variants, tps, fps, fns = check(vcfpy.Reader(args.merged_vcf), csv.DictReader(args.true_variants), csv.DictReader(args.variant_cov))
    fieldnames = ['CHROM','POS','REF','ALT','TYPE','INS_TYPE','LEN','QUAL','AD_REF','AD_ALT','DP','BAM_DP','GT','ZYG','RK_DISCOSNP','TOOL','SAMPLE','TRUE_POS','FALSE_POS','FALSE_NEG']
    write_variants(all_variants, fieldnames, args.called_variants)

if __name__ == '__main__':
    sys.exit(main())
