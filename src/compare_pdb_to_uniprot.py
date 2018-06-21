#!/usr/bin/env python3


from argparse import ArgumentParser
from Bio.PDB import PDBParser, is_aa, Polypeptide
import sys


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
    parser = ArgumentParser()
    parser.add_argument('-p', '--pdb', dest='pdb', required=True, help='PDB file')
    parser.add_argument('-f', '--fasta', dest='fasta', required=True, help='fasta file')
    parser.add_argument('-k', '--k', dest='k', type=int, default=10, help='check the first k residues in the PDB file')
    args = parser.parse_args()

    # read sequence
    sequence = read_fasta(args.fasta)

    # read PDB file
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(id='pdb', file=args.pdb)
    model = structure[0]
    residues = [r for r in model.get_residues() if is_aa(r, standard=True)]

    # check that PDB residue indices match UniProt residue indices 
    # for the firt 10 PDB residues
    pdb_seq = ''
    uniprot_seq = ''
    for r in residues[:args.k]:
        pdb_res_idx = r.get_id()[1]
        pdb_seq += Polypeptide.three_to_one(r.get_resname())
        uniprot_seq += sequence[pdb_res_idx-1]

    if pdb_seq == uniprot_seq:
        print('Residue indices in the PDB file match those in the '
              'UniProt fasta file for the first', args.k, 'residues.')
    else:
        print('Residue indices in the PDB file does not match those in the '
              'UniProt fasta file for the first', args.k, 'residues.')
        print('Check files', args.pdb, args.fasta)

    print('PDB    :', pdb_seq)
    print('UniProt:', uniprot_seq)


if __name__ == '__main__':
    main()
