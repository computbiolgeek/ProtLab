#!/usr/bin/env python3

import sys

def read_fasta(fasta_file):
    with open(fasta_file, 'rt') as f:
        fasta_seq = ''.join([l.strip() for l in f.readlines() if not l.strip().startswith(('>', '#'))])
    return fasta_seq


def one2three(input_code):
    one2three = {
        'A': 'ALA',
        'R': 'ARG',
        'N': 'ASN',
        'D': 'ASP',
        'C': 'CYS',
        'E': 'GLU',
        'Q': 'GLN',
        'G': 'GLY',
        'H': 'HIS',
        'I': 'ILE',
        'L': 'LEU',
        'K': 'LYS',
        'M': 'MET',
        'F': 'PHE',
        'P': 'PRO',
        'S': 'SER',
        'T': 'THR',
        'W': 'TRP',
        'Y': 'TYR',
        'V': 'VAL'
    }

    # make sure that input is valid
    assert len(input_code) == 1, 'input_code must have extactly one letter'
    assert input_code in one2three.keys(), 'given input_code was not recognized'

    # return the one-letter code
    input_code = input_code.upper()
    return one2three[input_code]


def main():
    '''
    '''
    args = sys.argv[1:]

    # print help info
    if '-h' in args:
        print('usage: fasta2seqres <fasta file> <chain id>')
        sys.exit(0)

    # parse the given fasta file
    fasta_seq = read_fasta(args[0])
    chain_id = args[1]

    # generate SEQRES lines
    seq_length = len(fasta_seq)
    res_per_line = 13
    seqres_lines = []
    for i in range(0, seq_length, res_per_line):
        if (i + res_per_line) <= seq_length:
            seqres_lines.append(' '.join([one2three(x) for x in fasta_seq[i:i+res_per_line]]))
        else:
            seqres_lines.append(' '.join([one2three(x) for x in fasta_seq[i:]]))
    
    #
    j = 1
    for line in seqres_lines:
        print('SEQRES  ' + '%2s' % str(j) + ' ' +  chain_id + ' ' + '%4s' % str(seq_length) + '  ' + line)
        j += 1



if __name__ == '__main__':
    main()
