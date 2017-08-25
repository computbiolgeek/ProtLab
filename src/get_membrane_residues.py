#!/usr/bin/env python

"""
    @summary: Given an OPM PDB file, this script extracts residues that are within the membrane,
    and returns them into a new PDB file.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
"""

from Bio.PDB import *
from argparse import ArgumentParser

# parse command line options and arguments
parser = ArgumentParser()
parser.add_argument("--input", help = "input pdb file downloaded from the OPM database")
parser.add_argument("--thickness", help = "1/2 of the membrane bilayer thickness, can usually \
    be found on the first line of the pdb file", default = 15, type = float)
parser.add_argument("--output", help = "output pdb file")
args = parser.parse_args()

class MemSelect(Select):
    """
        @summary: Select class which will be used to select transmembrane residues.
    """
    def __init__(self, thickness):
        """
            @summary: constructor
        """
        self.thickness = thickness
         
    def accept_residue(self, residue):
        """
            @summary: One of the four hook methods: accept_model(), accept_chain(), accept_residue(), and accept_atom()
        """
        if not 'CA' in residue:
            return False
        ca_coords = residue['CA'].get_coord()
        if ca_coords[2] > self.thickness or ca_coords[2] < -self.thickness:
            return False
        else:
            return True
    

# process the input pdb file and write the output
pdb_parser = PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure("input_pdb", args.input)
pdb_io = PDBIO()
pdb_io.set_structure(structure)
pdb_io.save(args.output, MemSelect(args.thickness))
