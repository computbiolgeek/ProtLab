#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
from Bio.PDB import PDBParser
import numpy as np
import pandas as pd

def compute_sc_centroid(residue):
    """
        C = (v1 + ... + vn) / n
    """
    # return CA coordinates if residue is GLY
    if residue.get_resname() == "GLY":
        return residue['CA'].get_coord()
    # for other residues, return sidechain centroid
    backbone_atoms = {'N', 'CA', 'C', 'O'}
    centroid = np.array( [0.0, 0.0, 0.0] )
    number_sidechain_atoms = 0
    for atom in residue.get_list():
        if atom.get_name() not in backbone_atoms:
            centroid += atom.get_coord()
            number_sidechain_atoms += 1
    return centroid / number_sidechain_atoms


def compute_distance(res_a, res_b):
    '''
    The Euclidean distance between the centroid of residue A and the centroid of residue B.
    '''
    centroid_a = compute_sc_centroid(res_a)
    centroid_b = compute_sc_centroid(res_b)
    return np.linalg.norm(centroid_a - centroid_b)


def main():
    '''
    '''
    # command line argument parser
    parser = ArgumentParser(description='Compute a distance matrix given two list of residues.')
    parser.add_argument('-p', '--pdb', dest='pdb', help='input file in PDB format')
    parser.add_argument('-r', '--res', dest='res', nargs=2, help='two files, each containing a list of residue IDs')
    parser.add_argument('-o', '--output', dest='output', help='name to the file to which to write the distance matrix')

    # parse command line arguments
    args = parser.parse_args()

    # parse the given PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='pdb', file=args.pdb)
    model = structure[0]

    # get files containing residues
    res_files = args.res

    # get the residues of the first chain 
    with open(res_files[0], 'rt') as f:
        res_ids_a = [tuple(line.strip().split()) for line in f.readlines()]
    chain_a = res_ids_a[0][0]
    residues_a = [model[chain_a][int(pos[1])] for pos in res_ids_a]

    # get the residues of the second chain
    with open(res_files[1], 'rt') as f:
        res_ids_b = [tuple(line.strip().split()) for line in f.readlines()]
    chain_b = res_ids_b[0][0]
    residues_b = [model[chain_b][int(pos[1])] for pos in res_ids_b]

    # compute pair distances between residues in the first chain and residues in the second chain
    pair_distances = np.array([[compute_distance(a, b) for b in residues_b] for a in residues_a])
    indices = [pos[1] for pos in res_ids_a]
    names = [pos[1] for pos in res_ids_b]
    pair_distances_df = pd.DataFrame(pair_distances, index=indices, columns=names)

    # write the matrix to file
    pair_distances_df.to_csv(args.output, index=True, header=True, sep=',', float_format='%.2f')


if __name__ == '__main__':
    main()

