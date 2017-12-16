#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
from Bio.PDB import PDBParser

def get_residues(model, chain_id, res_ids):
    '''
    '''
    # get the chain
    chain = model[chain_id]

    # return the requested residues
    for res_id in res_ids:
        yield chain[res_id]
        

def read_res_ids(res_file):
    '''
    '''
    #
    with open(res_file, 'rt') as f:
        # the first line is assumed to contain the chain ID
        chain_id = f.readline().strip()
        # the rest lines are assumed to be residue IDs
        res_ids = [int(line.strip()) for line in f.readlines()]
    return chain_id, res_ids



def main():
    '''
    '''
    # command line argument parser
    parser = ArgumentParser()
    parser.add_argument('--pdb', dest='pdb', help='input file in PDB format')
    parser.add_argument('--res', dest='res', help='two files, each containing a list of residue IDs')
    parser.add_argument('--output', dest='output', help='name to the file to which to write the distance matrix')

    # parse command line arguments
    args = parser.parse_args()

    # parse the given PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='pdb', file=args.pdb)
    model = structure[0]

