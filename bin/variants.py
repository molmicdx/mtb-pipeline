
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
                        help='reference genome to intoduce mutations into')
    parser.add_argument('--settings',
                        default="variants_settings.conf",
                        help="Settings file")
    parser.add_argument('--evendist', action='store_true', help='evenly distribute variants')
    parser.add_argument('mutations', type=argparse.FileType('w'),
                        help='file with list of introduced mutations')
    parser.add_argument('output', type=argparse.FileType('w'),
                        help='name of output file with mutated genome')
    return  parser.parse_args()


def add_ndist(mean, stddev, mindist):
    return max(int(abs(np.random.normal(mean, stddev))), mindist)

               
def add_dist(density):
    return randint(1, 1/density*2+1)


def main():
    args = get_args()
    settings = ConfigParser(allow_no_value=False)
    settings.read(args.settings)
    reference  = next(fastalite(args.reference))
    newseq = list(reference.seq)
    name = reference.id

    relpos = -1
    abspos = -1
    mutations = []
    snps = 0
    indels = []

    while True:
        # even distribution of variants
        if(args.evendist):
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

            mutations.append([name, abspos + 1, '.', old, snp])
            snps += 1
        else:
            indel_len = randint(settings.getint('variants', 'indelmin'),
                                settings.getint('variants', 'indelmax'))
            # choose an insertion vs a deletion
            if random() < settings.getfloat('variants', 'insertion_to_deletion'):
                char_at = newseq[relpos]
                # choose indel from forward or back based on availability
                insertion = ''
                if(relpos + 1 + indel_len > len(newseq) - 1):
                    continue

                    # TODO: remove this?
#                insordel = -1
#                if(relpos + indel_len > len(newseq)):
#                    insordel = 1
#                elif(relpos - indel_len < 0):
#                    insordel = 0
#                else:
#                    insordel = random()

                #                if insordel > settings.getfloat('variants', 'insertion_to_deletion'):
                insertion = newseq[relpos+1:relpos+1+indel_len]
#                else:
#                    insertion = newseq[relpos-indel_len:relpos]

                newseq[relpos+1:relpos+1] = insertion

                mutations.append([name, abspos + 1, '.', char_at, char_at + ''.join(insertion)])
                relpos += indel_len
                indels.append(('insertion', len(insertion)))
            else:
                char_at = newseq[relpos]
                removed = newseq[relpos+1:relpos + 1 + indel_len]
                newseq[relpos+1:relpos + 1 + indel_len] = []
                mutations.append([name, abspos + 1, '.', char_at + ''.join(removed), char_at])
                abspos += indel_len
                indels.append(('deletion', len(removed)))

    args.output.write('>' + reference.id + '\n')
    args.output.write(''.join(newseq))

    writer = csv.writer(args.mutations, delimiter='\t')
    writer.writerow(['#CHROM', 'POS', 'ID', 'REF', 'ALT'])
    writer.writerows(mutations)

    print('tempate length:', len(reference.seq))
    print('output length:', len(newseq))
    print('number of mutations', len(mutations))
    print('number of snps:', snps)
    # number and length distribution of insertions and deletions (create pandas data.frame from 'stats'?)

if __name__ == '__main__':
    sys.exit(main())
