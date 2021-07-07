import argparse
import csv
import sys
import vcfpy


def get_args():
    parser = argparse.ArgumentParser(description="Check pipeline results with true variants. Outputs 3 CSVs: false positives, false negatives, and summary including precision, recall")
    parser.add_argument('merged_vcf', type=argparse.FileType('r'),
                        help="merged vcf file with pipeline results")
    parser.add_argument('true_variants', type=argparse.FileType('r'),
                        help="csv file with sorted true variants")
    parser.add_argument('false_positives', type=argparse.FileType('w'),
                        help="output csv with false posistive variants")
    parser.add_argument('false_negatives', type=argparse.FileType('w'),
                        help="output csv with false negative variants")
    return parser.parse_args()


def get_variant_from_vcf_record(record):
    variant = {}
    variant['CHROM'] = record.CHROM
    variant['POS'] = str(record.POS)
    variant['REF'] = record.REF
    variant['ALT'] = ','.join([alt.value for alt in record.ALT])
    for call in record.calls:
        variant[call.sample] = call.data['GT']
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
    alt = vcf_record.ALT == true_variant['ALT']
    #type = ','.join(vcf_record.INFO['TY']) == true_variant['TYPE']
    genotypes = True
    for call in vcf_record.calls:
        try: #Catch diploid GT reporting
            if call.data['GT'][0] != true_variant[call.sample] and call.data['GT'][2] != true_variant[call.sample]:
                genotypes = False
        except IndexError:
            if call.data['GT'][0] and call.data['GT'][0] != true_variant[call.sample]:
                genotypes = False
    return chrom and pos and ref and alt and genotypes


def check(vcf_reader, true_variants_reader):
    tps, fps, fns = [], [], []
    vcf_record = next(vcf_reader, None)
    true_variant = next(true_variants_reader, None)
    while vcf_record and true_variant:
        if (vcf_record.CHROM, vcf_record.POS) <= (true_variant['CHROM'], int(true_variant['POS'])):
            if is_match(vcf_record, true_variant):
                tps.append(true_variant)
                true_variant = next(true_variants_reader, None)
            else:
                fps.append(get_variant_from_vcf_record(vcf_record))
            vcf_record = next(vcf_reader, None)
        else:
            fns.append(true_variant)
            true_variant = next(true_variants_reader, None)
    while vcf_record:
        fps.append(get_variant_from_vcf_record(vcf_record))
        vcf_record = next(vcf_reader, None)
    while true_variant:
        fns.append(true_variant)
        true_variant = next(true_variants_reader, None)
    return tps, fps, fns

def stats(tps, fps, fns):
    tp = str(len(tps))
    fp = str(len(fps))
    fn = str(len(fns))
    
    snp = [0,0,0]
    ins = [0,0,0]
    dele = [0,0,0]
    for variant in tps:
        if variant['TYPE'] == 'SNP':
            snp[0] += 1
        elif variant['TYPE'] == 'INS':
            ins[0] += 1
        elif variant['TYPE'] == 'DEL':
            dele[0] += 1
    for variant in fps:
        if variant['TYPE'] == 'SNP':
            snp[1] += 1
        elif variant['TYPE'] == 'INS':
            ins[1] += 1
        elif variant['TYPE'] == 'DEL':
            dele[1] += 1
    for variant in fns:
        if variant['TYPE'] == 'SNP':
            snp[2] += 1
        elif variant['TYPE'] == 'INS':
            snp[2] += 1
        elif variant['TYPE'] == 'DEL':
            snp[2] += 1
    try:
        precision = str(len(tps)/(len(tps) + len(fps)))
    except ZeroDivisionError:
        precision = '.'
    try:
        recall = str(len(tps)/(len(tps) + len(fns)))
    except ZeroDivisionError:
        recall = '.'

    return tp, fp, fn, snp, ins, dele, precision, recall

def write_variants(variants, fieldnames, file):
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(variants)


def main():
    args = get_args()
    tps, fps, fns = check(vcfpy.Reader(args.merged_vcf), csv.DictReader(args.true_variants))
    try: #catch exception when using reference reads, i.e., no true variants
        fieldnames = tps[0].keys()
    except IndexError:
        with open(args.true_variants.name, 'r') as tvfile:
            header = tvfile.readline().rstrip()
        fieldnames = {col: ',' for col in header.split(',')}.keys()
    write_variants(fps, fieldnames, args.false_positives)
    write_variants(fns, fieldnames, args.false_negatives)
    tp, fp, fn, snp, ins, dele, precision, recall = stats(tps, fps, fns)
    output_dir = '/'.join(args.false_positives.name.split('/')[:-1])
    samplename = args.merged_vcf.name.split('/')[-1].split('.')[0]
    with open(output_dir + '/' + samplename + '_stats.csv', 'w') as vcfile:
        vcfile.write('SAMPLE,TP,TP_SNP,TP_IND,FP,FP_SNP,FP_IND,FN,FN_SNP,FN_IND,PRECISION,RECALL\n')
        vcfile.write(samplename + ',' \
                     + tp + ',' + str(snp[0]) + ',' + str(ins[0] + dele[0]) + ',' \
                     + fp + ',' + str(snp[1]) + ',' + str(ins[1] + dele[1]) + ',' \
                     + fn + ',' + str(snp[2]) + ',' + str(ins[2] + dele[2]) + ',' \
                     + precision + ',' + recall + '\n')

if __name__ == '__main__':
    sys.exit(main())
