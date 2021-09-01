import argparse
import csv
import sys

parser = argparse.ArgumentParser(description='Add read depth to variants csv file. Input bed file is produced by bedtools from the intersection of variants bed and genome coverage bed files, with number of bp overlap in the final column')

parser.add_argument('variants_bed', type=argparse.FileType('r'), help='variant bed file with genome coverage information')
parser.add_argument('variants_csv', type=argparse.FileType('r'), help='csv file of variants')
parser.add_argument('variants_cov', type=argparse.FileType('w'), help='output csv with variants and average read depths (i.e. BAM_DP)')
parser.add_argument('--bam_dp_out', type=argparse.FileType('w'), default=None, help='output csv with positions, entry or mutation type, and average read depths (i.e. BAM_DP)')

def get_read_depths(var_bed_in):
    var_reader = csv.reader(var_bed_in, delimiter='\t')
    mutation_bed = next(var_reader, None)
    cov_along_mutation = {}
    while mutation_bed:
        # cov dictionary key is (chromStart + 1): value is list of tuples of (read depth, number of bps, left flank chromStart, right flank chromStart).
        # BED uses zero-based positions, VCF (and variants.py generated files) uses 1-based positions
        csv_POS = int(mutation_bed[1]) + 1
        read_depth = int(mutation_bed[-2])
        num_bp = int(mutation_bed[-1])
        entry_type = mutation_bed[3]
        lflankpos = mutation_bed[4]
        rflankpos = mutation_bed[5]
        try:
            cov_along_mutation[csv_POS].append((read_depth, num_bp, lflankpos, rflankpos, entry_type))
        except KeyError:
            cov_along_mutation[csv_POS] = [(read_depth, num_bp, lflankpos, rflankpos, entry_type)]
        mutation_bed = next(csv.reader(var_bed_in, delimiter='\t'), None)
    return cov_along_mutation

def calculate_bam_dp(read_depth_list):
    read_depths_sum = 0
    total_bp = 0
    for entry in read_depth_list:
        read_depths_sum += entry[0] * entry[1]
        total_bp += entry[1]
    bam_dp = round(read_depths_sum / total_bp)
    return bam_dp

def get_variants_bam_dp(cov_along_mutation, var_bed_in):
    variants_bam_dp = {}
    for pos in cov_along_mutation.keys():
        counter = 0
        while counter < len(cov_along_mutation[pos]):
            if (cov_along_mutation[pos][counter][2] == '.') or (cov_along_mutation[pos][counter][3] == '.'):
                variants_bam_dp[pos] = [calculate_bam_dp(cov_along_mutation[pos]), cov_along_mutation[pos][counter][-1]]
                counter += 1
            else:
                lflankpos = cov_along_mutation[pos][counter][2]
                rflankpos = cov_along_mutation[pos][counter][3]
                left_bam_dp = calculate_bam_dp(cov_along_mutation[int(lflankpos) + 1])
                right_bam_dp = calculate_bam_dp(cov_along_mutation[int(rflankpos) + 1])
                variants_bam_dp[pos] = [round((left_bam_dp + right_bam_dp)/2), cov_along_mutation[pos][counter][-1]]
                # troubleshooting lancet edge case where lflankpos=1047759 is also chromStart for another DEL
                # resulting in extraneous values provided to calculate_bam_dp(). Not fixed yet.
                #if pos == 1047765:
                    #print(cov_along_mutation[int(lflankpos) + 1])
                    #print(lflankpos, rflankpos, left_bam_dp, right_bam_dp)
                break
    return variants_bam_dp

def write_cov_to_csv(var_bam_dp, var_csv_in, cov_csv_out):
    variants_reader = csv.DictReader(var_csv_in)
    writer = csv.DictWriter(cov_csv_out, fieldnames = variants_reader.fieldnames)
    writer.writeheader()
    variant = next(variants_reader, None)
    while variant:
        variant['BAM_DP'] = var_bam_dp[int(variant['POS'])][0]
        writer.writerow(variant)
        variant = next(variants_reader, None)
    return

def get_original_bam_dp(cov_along_mutation, var_bed_in):
    original_bam_dp = {}
    for pos in cov_along_mutation.keys():
        counter = 0
        while counter < len(cov_along_mutation[pos]):
            original_bam_dp[pos] = [calculate_bam_dp(cov_along_mutation[pos]), cov_along_mutation[pos][counter][-1]]
            counter += 1
    return original_bam_dp

def write_original_bam_dp(orig_bam_dp, bamdp_csv_out):
    entry = {}
    fieldnames = ['POS','TYPE','BAM_DP']
    writer = csv.DictWriter(bamdp_csv_out, fieldnames=fieldnames)
    writer.writeheader()
    for pos in orig_bam_dp.keys():
        entry['POS'] = pos
        entry['TYPE'] = orig_bam_dp[pos][1]
        entry['BAM_DP'] = orig_bam_dp[pos][0]
        writer.writerow(entry)
    return

def main():
    args = parser.parse_args()
    cov_at_mutation_pos = get_read_depths(args.variants_bed)
    mutation_bam_dp = get_variants_bam_dp(cov_at_mutation_pos, args.variants_bed)
    write_cov_to_csv(mutation_bam_dp, args.variants_csv, args.variants_cov)
    if args.bam_dp_out != None:
        original_bam_dp = get_original_bam_dp(cov_at_mutation_pos, args.variants_bed)
        write_original_bam_dp(original_bam_dp, args.bam_dp_out)

if __name__ == '__main__':
    sys.exit(main())

