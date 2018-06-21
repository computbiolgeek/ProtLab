#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import os, re


def read_fasta(fastafile):
    '''
    Parse a fasta file into a list of dicts, each dict is a 
    '''
    seq = ''
    with open(fastafile, 'rt') as f:
        # get the first line
        first_line = next(f)
        if not first_line.startswith('>'):
            print('Invalid fasta file, first line must start with \'>\'.')
            sys.exit(1)

        # read the remaining lines
        for line in f:
            # skip empty lines
            if not line.strip():
                continue
            if line.startswith('>'):
                # found a new sequence record
                print('More than one record found, please remove extra records.')
                sys.exit(1)
            else:
                # concatenate sequence lines
                seq += line.strip()
    return seq


def main():
    '''
    '''
    # parse command line arguments
    parser = ArgumentParser()
    parser.add_argument('-f', '--fasta', dest='fasta', help='input FASTA file')
    parser.add_argument('-r', '--ranges', dest='ranges', help='transmembrane helix residue number ranges')
    parser.add_argument('-o', '--output', dest='output', help='output FASTA file')
    args = parser.parse_args()

    # extract sequences from residues with resolved coordinates
    sequence = read_fasta(args.fasta)
    fasta_id = os.path.basename(args.fasta).split('.')[0]

    # parse the transmembrane helix residue number ranges
    with open(args.ranges, 'rt') as f:
        tmh_ranges = [tuple(re.split(r'[:-]', l.strip())) for l in f]

    # write extracted transmembrane helix sequences into a fasta file
    tmh_seqs = ['> ' + fasta_id + ' TMH sequences\n']
    with open(args.output, 'wt') as f:
        for r in tmh_ranges:
            start_index = int(r[1]) - 1
            end_index = int(r[2]) - 1
            tmh_seqs.append(sequence[start_index:end_index+1] + '\n')
        f.writelines(tmh_seqs)
        

if __name__ == '__main__':
    main()
