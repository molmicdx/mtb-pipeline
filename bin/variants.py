
#!/usr/bin/env python3

"""Simulate mutations in given genome.

"""

import sys
from random import randint, choice, random
import argparse
import collections
import csv
from configparser import ConfigParser

import numpy as np
from fastalite import fastalite, Opener

ALPHABET="ACGT"


def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Opener('r'),
                        help='reference genome to introduce mutations into')
    parser.add_argument('--settings',
                        default="configs/variants_settings.conf",
                        help="Settings file")
    parser.add_argument('--evendist', action='store_true', help='evenly distribute variants')
    parser.add_argument('--inputmutations', type=argparse.FileType('r'), default=None, help='csv of mutations to introduce')
    parser.add_argument('mutations', type=argparse.FileType('w'),
                        help='file with list of introduced mutations')
    parser.add_argument('output', type=argparse.FileType('w'),
                        help='name of output file with mutated genome')
    return  parser.parse_args()


def add_ndist(mean, stddev, mindist):
    return max(int(abs(np.random.normal(mean, stddev))), mindist)

               
def add_dist(density):
    return randint(1, 1/density*2+1)


def add_specific_mutations(mutations_csv, ref_genome, mutations_out, var_genome_fa):
    mutations_reader = csv.DictReader(mutations_csv)
    reference = next(fastalite(ref_genome))
    ref_seq = list(reference.seq)
    print(len(ref_seq))
    new_seq = []
    name = reference.id
    fieldnames = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'TYPE', 'INS_TYPE']
    writer = csv.DictWriter(mutations_out, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    mutation = next(mutations_reader, None)
    start_pos = 0
    while mutation:
        end_pos = int(mutation['POS'])
        print(start_pos, end_pos)
        new_seq += ref_seq[start_pos:end_pos - 1]
        new_mutation = {}
        new_mutation['#CHROM'] = name
        new_mutation['POS'] = mutation['POS'] # this may change after identification of core genome set.
        new_mutation['ID'] = '.'
        new_mutation['TYPE'] = mutation['TYPE']
        new_mutation['INS_TYPE'] = mutation['INS_TYPE']
        #intro SNP into new_seq
        if new_mutation['TYPE'] == 'SNP':
            new_mutation['REF'] = ref_seq[end_pos - 1]
            if new_mutation['REF'] == mutation['REF']:
                new_mutation['ALT'] = mutation['ALT']
            else:
                new_mutation['ALT'] = choice(ALPHABET.replace(ref_seq[end_pos - 1], ''))
            start_pos = end_pos
        #intro DEL into newseq
        elif new_mutation['TYPE'] == 'DEL':
            new_mutation['REF'] = ''.join(ref_seq[end_pos - 1:(end_pos - 1 + len(mutation['REF']))])
            new_mutation['ALT'] = ref_seq[end_pos - 1]
            start_pos = end_pos - 1 + len(mutation['REF'])
        #intro INS into newseq
        elif new_mutation['TYPE'] == 'INS':
            new_mutation['REF'] = ref_seq[end_pos - 1]
            if new_mutation['INS_TYPE'] == 'RDM':
                new_mutation['ALT'] = mutation['ALT']
            elif new_mutation['INS_TYPE'] == 'DUP':
                new_mutation['ALT'] = ''.join(ref_seq[end_pos - 1:(end_pos - 1 + len(mutation['ALT']))])
            elif new_mutation['INS_TYPE'] == 'INV':
                complemented = []
                for nt in ref_seq[end_pos:(end_pos + len(mutation['ALT']) - 1)]:
                    if nt == 'A':
                        complemented.append('T')
                    elif nt == 'C':
                        complemented.append('G')
                    elif nt == 'T':
                        complemented.append('A')
                    elif nt == 'G':
                        complemented.append('C')
                new_mutation['ALT'] = ref_seq[end_pos - 1] + ''.join(reversed(complemented))
            start_pos = end_pos
        new_seq.append(new_mutation['ALT'])
        writer.writerow(new_mutation)
        mutation = next(mutations_reader, None)
    new_seq += ref_seq[start_pos:len(ref_seq)]
    var_genome_fa.write('>'+ name + ' variant\n')
    var_genome_fa.write(''.join(new_seq))
    print(len(''.join(new_seq)))
    return


def simulate_mutations(ref_genome, settings_file, mutations_out, var_genome_fa, evendist=True):
    settings = ConfigParser(allow_no_value=False)
    settings.read(settings_file)
    reference = next(fastalite(ref_genome))
    newseq = list(reference.seq)
    name = reference.id
    relpos = -1
    abspos = -1
    mutations = []
    snps = 0
    indels = []
    num_dup = 0
    num_rdm = 0
    num_inv = 0
    num_del = 0
    while True:
        # even distribution of variants
        if evendist:
            dist = add_dist(settings.get('variants', 'density'))
        # normal distribution of variants
        else:
            dist = add_ndist(1/settings.getfloat('variants', 'density'),
                             settings.getint('normal_dist', 'stddev'),
                             settings.getint('normal_dist', 'mindist'))
        relpos += dist
        abspos += dist
        if(relpos >= len(newseq)):
            break

        # choose a SNP vs an indel
        if random() < settings.getfloat('variants', 'snptoindel'):
            old = newseq[relpos]
            snp = choice(ALPHABET.replace(old, ''))
            newseq[relpos] = snp

            mutations.append([name, abspos + 1, '.', old, snp, 'SNP'])
            snps += 1
        else:
            indel_len = randint(settings.getint('variants', 'indelmin'),
                                settings.getint('variants', 'indelmax'))
            # choose an insertion vs a deletion
            if random() < settings.getfloat('variants', 'insertion_to_deletion'):
                char_at = newseq[relpos]
                # choose indel from forward or back based on availability
                if(relpos + 1 + indel_len > len(newseq) - 1):
                        continue
                # duplications
                if random() < settings.getfloat('variants', 'duplication_to_other_ins'):
                    insertion = newseq[relpos+1:relpos+1+indel_len]
                    mutations.append([name, abspos + 1, '.', char_at, char_at + ''.join(insertion), 'INS', 'DUP'])
                    indels.append(('duplication', len(insertion)))
                    num_dup += 1
                else:
                    insertion = []
                    # equal probability of random insertion sequence or inversion
                    if random() < 0.5:
                        while len(insertion) < indel_len:
                            insertion.append(ALPHABET[randint(0,3)])
                        mutations.append([name, abspos + 1, '.', char_at, char_at + ''.join(insertion), 'INS', 'RDM'])
                        indels.append(('random_ins', len(insertion)))
                        num_rdm += 1
                    else:
                        for nt in reversed(newseq[relpos+1:relpos+1+indel_len]):
                            if nt == 'A':
                                insertion.append('T')
                            elif nt == 'T':
                                insertion.append('A')
                            elif nt == 'C':
                                insertion.append('G')
                            elif nt == 'G':
                                insertion.append('C')
                        mutations.append([name, abspos + 1, '.', char_at, char_at + ''.join(insertion), 'INS', 'INV'])
                        indels.append(('inversion', len(insertion)))
                        num_inv += 1

                newseq[relpos+1:relpos+1] = insertion
                relpos += indel_len
            else:
                char_at = newseq[relpos]
                removed = newseq[relpos+1:relpos + 1 + indel_len]
                newseq[relpos+1:relpos + 1 + indel_len] = []
                mutations.append([name, abspos + 1, '.', char_at + ''.join(removed), char_at, 'DEL'])
                abspos += indel_len
                indels.append(('deletion', len(removed)))
                num_del += 1
    var_genome_fa.write('>' + reference.id + '\n')
    var_genome_fa.write(''.join(newseq))

    writer = csv.writer(mutations_out, delimiter='\t')
    writer.writerow(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'TYPE', 'INS_TYPE'])
    writer.writerows(mutations)
    
    print('Template length:', len(reference.seq))
    print('Output genome length:', len(newseq))
    print('Number of mutations (SNP + indels):', len(mutations))
    print('Number of SNP:', snps)
    print('Number of indels:', len(indels))
    print('Number of duplications (DUP):', num_dup)
    print('Number of random-sequence insertions (RDM):', num_rdm)
    print('Number of inversions (INV):', num_inv)
    print('Number of deletions (DEL):', num_del)
    return

def main():
    args = get_args()
    if args.inputmutations == None:
        simulate_mutations(args.reference, args.settings, args.mutations, args.output, evendist=args.evendist)
    else:
        add_specific_mutations(args.inputmutations, args.reference, args.mutations, args.output)
    return

if __name__ == '__main__':
    sys.exit(main())
