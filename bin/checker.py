import argparse
import csv
import sys
import vcfpy


def get_args():
    parser = argparse.ArgumentParser(description="Check pipeline results with true variants")
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
    variant = {call.sample: call.data['GT'] for call in record.calls}
    variant['CHROM'] = record.CHROM
    variant['POS'] = str(record.POS)
    variant['REF'] = record.REF
    variant['ALT'] = ','.join([alt.value for alt in record.ALT])
    variant['TYPE'] = ','.join(record.INFO['TY'])
    return variant


def is_match(vcf_record, true_variant):
    chrom = vcf_record.CHROM == true_variant['CHROM']
    pos = vcf_record.POS == int(true_variant['POS'])
    type = ','.join(vcf_record.INFO['TY']) == true_variant['TYPE']
    genotypes = True
    for call in vcf_record.calls:
        if call.data['GT'] and call.data['GT'] != true_variant[call.sample]:
            genotypes = False
    return chrom and pos and type and genotypes


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


def write_variants(variants, fieldnames, file):
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(variants)


def main():
    args = get_args()
    tps, fps, fns = check(vcfpy.Reader(args.merged_vcf), csv.DictReader(args.true_variants))
    fieldnames = tps[0].keys()
    write_variants(fps, fieldnames, args.false_positives)
    write_variants(fns, fieldnames, args.false_negatives)


if __name__ == '__main__':
    sys.exit(main())



