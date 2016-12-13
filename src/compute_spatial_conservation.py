"""
    @summary: 
    @author: 
    @contact: 
    @copyright: 
    @change: 
"""

from argparse import ArgumentParser
from Bio.PDB import PDBParser, NeighborSearch

parser = ArgumentParser( description = "Hi there, I compute residue spatial conservation for you!" )
parser.add_argument( "-p", "--pdb", dest = "pdb", required = True,
                    help = "PDB file" )
parser.add_argument( "-c", "--conservation", dest = "conservation", required = True,
                     help = "file that stores residue conservations computed from multiple sequence alignment" )
parser.add_argument( "-o", "outfile", dest = "outfile", required = True,
                    help = "output file name" )
parser.add_argument( "-v", "--verbose", dest = "verbose", action = "store_true",
                    help = "verbose mode" )
args = parser.parse_args()




