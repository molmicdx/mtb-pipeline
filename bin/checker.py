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
    #mutation_types = ['TYPE', 'Ty'] 
    #for mut in mutation_types:
    #    if mut in record.INFO.keys():
    #        variant['TYPE'] = ','.join(record.INFO[mut])
    if variant['REF'].startswith('<') or variant['ALT'].startswith('<'):
        variant['TYPE'] = variant['ALT'][1:-1] # delly reports type in ALT col
    else:
        if len(variant['REF']) == 1 and len(variant['ALT']) == 1:
            variant['TYPE'] = 'SNP'
        elif len(variant['REF']) > len(variant['ALT']):
            variant['TYPE'] = 'DEL'
        elif len(variant['REF']) < len(variant['ALT']):
            variant['TYPE'] = 'INS'
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


def check(vcf_reader, true_variants_reader):
    all_variants = []
    tps, fps, fns = [], [], []
    vcf_record = next(vcf_reader, None)
    true_variant = next(true_variants_reader, None)
    while vcf_record and true_variant:
        called = get_variant_from_vcf_record(vcf_record)
        if (vcf_record.CHROM, vcf_record.POS) <= (true_variant['CHROM'], int(true_variant['POS'])):
            if is_match(vcf_record, true_variant):
                called['TRUE_POS'] = 1
                all_variants.append(called)
                tps.append(called)
                true_variant = next(true_variants_reader, None)
            else:
                called['FALSE_POS'] = 1
                all_variants.append(called)
                fps.append(called)
            vcf_record = next(vcf_reader, None)
        else:
            all_variants.append(true_variant)
            fns.append(true_variant)
            true_variant = next(true_variants_reader, None)
    while vcf_record:
        called['FALSE_POS'] = 1
        all_variants.append(called)
        fps.append(called)
        vcf_record = next(vcf_reader, None)
    while true_variant:
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
    all_variants, tps, fps, fns = check(vcfpy.Reader(args.merged_vcf), csv.DictReader(args.true_variants))
    fieldnames = ['CHROM','POS','REF','ALT','TYPE','QUAL','AD_REF','AD_ALT','DP','BAM_DP','GT','RK_DISCOSNP','TRUE_POS','FALSE_POS','FALSE_NEG','TOOL','SAMPLE']
    write_variants(all_variants, fieldnames, args.called_variants)

if __name__ == '__main__':
    sys.exit(main())
