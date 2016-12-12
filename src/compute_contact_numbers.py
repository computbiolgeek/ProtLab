"""
    @summary:
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @copyright: Bian Li
    @change:
"""

from argparse import ArgumentParser
from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import csv

parser = ArgumentParser( description = "Hi there, I compute contact numbers for you." )
parser.add_argument( "-i", "--infile", dest = "infile", required = True,
                    help = "input PDB file for computing contact numbers" )
parser.add_argument( "-o", "--outfile", dest = "outfile", required = True,
                    help = "output file to store computed contact numbers" )
parser.add_argument( "-t", "--type", dest = "type", choices = {"weighted", "unweighted"},
                     default = "weighted", help = "the type of contact number definition" )
parser.add_argument( "-b", "--bounds", dest = "bounds", nargs = 2, type = float, required = False,
                     default = [4, 11.4], help = "cutoff distances within which residues will be considered" )
parser.add_argument( "-s", "--sequence_separation", dest = "sequence_separation", type = int,
                    required = False, default = 3, help = "sequence separation beyond which "
                    "residues will be considered" )
parser.add_argument( "-v", "--verbose", dest = "verbose", required = False, action = "store_true",
                    help = "verbose mode" )
args = parser.parse_args()

# print collected command line arguments
if args.verbose:
    print( "input file:          " + args.infile )
    print( "output file:         " + args.outfile )
    print( "type:                " + args.type )
    print( "bounds:              " + str( args.bounds ) )
    print( "sequence separation: " + str( args.sequence_separation ) )
    print( "verbose:             " + str( args.verbose ) )

# parse the given PDB file
pdb_parser = PDBParser( PERMISSIVE = 1 )
structure = pdb_parser.get_structure( "pdb", args.infile )

def ComputeSideChainCentroid( residue ):
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
    for atom in residue.get_atom():
        if atom.get_name() not in backbone_atoms:
            centroid += atom.get_coord()
            number_sidechain_atoms += 1
    return centroid / number_sidechain_atoms

def ComputeContactNumber( residue, neighbor_list, bounds, cn_type = "unweighted" ):
    """
    """
    # compute contact number
    cn = 0
    if cn_type == "unweighted":
        cn = len( neighbor_list )
    else:
        centroid = ComputeSideChainCentroid( residue )
        for neighbor in neighbor_list:
            neighor_centroid = ComputeSideChainCentroid( neighbor )
            d_ij = np.linalg.norm( centroid - neighor_centroid )
            if d_ij <= bounds[0]:
                cn += 1
            else:
                cn += ( 1.0 / 2 ) * np.cos( ( d_ij - bounds[0] ) / ( bounds[1] - bounds[0] ) ) + 1.0 / 2
    return cn

# compute contact number for each residue
atom_list = list( structure.get_atoms() )
ns = NeighborSearch( atom_list )
contact_numbers = []
for residue in structure.get_residues():
    centroid = ComputeSideChainCentroid( residue )
    # compile a list of neighbors for the current residue
    neighbor_list = [res for res in ns.search( centroid, args.bounds[1], 'R' )
                     if res.get_id()[1] - residue.get_id()[1] > args.sequence_separation
                        or res.get_id()[1] - residue.get_id()[1] < -args.sequence_separation]
    residue_id = residue.get_id()[1]
    contact_numbers.append( 
        ( residue_id, ComputeContactNumber( residue, neighbor_list, args.bounds, args.type ) )
    )

# write contact numbers into a .csv file
with open( args.outfile, "wt" ) as f:
    csv_f = csv.writer( f )
    for row in contact_numbers:
        csv_f.writerow( row )
