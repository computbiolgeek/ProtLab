"""
    @summary: 
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @copyright: Bian Li
    @change: 
"""

from argparse import ArgumentParser
from Bio.PDB import PDBParser, NeighborSearch

parser = ArgumentParser()
parser.add_argument( "-i", "--infile", dest = "infile", required = True,
                    help = "input PDB file for computing contact numbers" )
parser.add_argument( "-o", "outfile", dest = "outfile", required = True,
                    help = "output file to store computed contact numbers" )
parser.add_argument( "-t", "--type", dest = "type", choices = {"weighted", "unweighted"},
                     default = "weighted", help = "the type of contact number definition" )
parser.add_argument( "-c", "--cutoff", dest = "cutoff", type = float, required = False,
                     default = 5, help = "cutoff distance within which residues will be considered" )
parser.add_argument( "-s", "--sequence_separation", dest = "sequence_separation", type = int,
                    required = False, default = 3, help = "sequence separation beyond which "
                    "residues will be considered" )
args = parser.parse_args()
