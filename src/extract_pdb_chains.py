#!/usr/bin/env python3

import sys


def main():
    '''
    '''
    args = sys.argv

    # if -h is given on the command line, print help information
    if '-h' in args:
        print('usage: extract_pdb_chains.py -c <chain IDs> -i <input PDB file> -o <output PDB file>')
        sys.exit(0)

    # read in chain IDs
    chain_ids = args[2]

    # read in the input PDB file
    with open(args[4], 'rt') as f:
        pdb_lines = [l for l in f.readlines() if l.startswith('ATOM') and l[21] in chain_ids]

    # write to output PDB file
    with open(args[6], 'wt') as f:
        f.write(''.join(pdb_lines))


if __name__ == '__main__':
    main()
