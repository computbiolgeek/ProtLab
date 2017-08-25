#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 09:48:34 2017

@author: lib14
"""

from argparse import ArgumentParser
from Bio.PDB import PDBParser, PDBIO


def main():
    """
    """
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True, help='input PDB file')
    parser.add_argument('-r', '--ranges', dest='ranges', required=True, 
                        help='ranges to be used to re-number residues')
    parser.add_argument('-c', '--chains', dest='chains', required=True,
                        help='chains to be considered')
    parser.add_argument('-o', '--output', dest='output', required=True, help='output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode')
    
    args = parser.parse_args()
    
    # print collected command line arguments
    if args.verbose:
        print('input file:  ' + args.input)
        print('output file: ' + args.output)
        print('range:       ' + str(args.range))
        print('start:       ' + str(args.start))
        print('verbose:     ' + str(args.verbose))
    
    # parse the given PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure('pdb', args.input)

    #
    ranges = args.ranges.split(',')
    seq = []
    for r in ranges:
        s, e = r.split('-')
        seq += list(range(int(s), int(e) + 1))
    
    for c in args.chains:
        i = 0
        for residue in structure[0][c]:
            residue.id = (' ', seq[i], ' ')
            i += 1
            if i == len(seq):
                break
            
    # write a new PDB file
    w = PDBIO()
    w.set_structure(structure)
    w.save(args.output)
            
if __name__ == '__main__':
    main()