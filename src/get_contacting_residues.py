#!/usr/bin/python

"""
    @summary: This script takes in a pdb file and a residue id. It computes the distance between
    the beta carbon atom of the residue of interest and that of all other residues. It finaly writes
    a list of residue ids along with the distances for residues within a specified cutoff.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
"""

from Bio.PDB import *
from argparse import ArgumentParser

# parse command line options and arguments
parser = ArgumentParser()
parser.add_argument("--pdb", help = "protein coordinate file in the PDB format")
parser.add_argument("--resi", help = "residue id", type = int)
parser.add_argument("--cutoff", help = "cutoff within which residues are considered", default = 10, type = float)
parser.add_argument("--output", help = "name of the output file", default = "contacting_residues.out")
args = parser.parse_args()

# parse in the pdb file
pdb_parser = PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure("pdb", args.pdb)

# get the beta carbon of interest
target_cb = structure[0]['A'][args.resi]['CB']
for residue in structure.get_residues():
    if not is_aa(residue, True) or (residue.get_id()[1] >= args.resi - 5 
                                    and residue.get_id()[1] <= args.resi + 5):
        continue
    else:
        if residue.has_id('CB'):
            cb_cb_dist = target_cb - residue['CB']
        else:
            cb_cb_dist = target_cb - residue['CA']
    if cb_cb_dist <= args.cutoff:
        resi_dist = (residue.get_id()[1], cb_cb_dist)
        print resi_dist