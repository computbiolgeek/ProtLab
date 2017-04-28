"""
    @summary:
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu or comput.biol.geek@gmail.com
    @change:
    @copyright:
"""

from argparse import ArgumentParser

from Bio.PDB import PDBParser

import matplotlib.pyplot as plt
import numpy as np


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

def ComputeCentroidDistance( res_a, res_b ):
    """
        @summary: Computes the distance between two residues.

        @param res_a: Biopython Residue object
        @param res_b: Biopython Residue object

        @return: Euclidean distance between res_a and res_b
    """
    centroid_a = ComputeSideChainCentroid( res_a )
    centroid_b = ComputeSideChainCentroid( res_b )
    return np.linalg.norm( centroid_a - centroid_b )

def ComputeRadiusOfGyration( coordinates ):
    """
        @summary: Computes the radius of gyration of the passed set of atoms.
        @param
            coordinates: a (N, 3) dimensional NumPy array.
        @return:
            Rg: a float with the corresponding radius of gyration
    """
    r_cm = coordinates.mean( axis = 0 )
    rg = np.mean( np.sum( np.linalg.norm( coordinates - r_cm, ord = 2, axis = 1 ) ) )
    return np.sqrt( rg )

def ComputeNormalizationFunction( r, a, r_c ):
    """
        @see: Shen and Sali, Protein Science, 2006
    """
    if r_c > 2 * a:
        return 3 * ( r * ( 2 - 2 * a ) ) ** 2 * ( r + 4 * a ) / ( 16 * a ** 6 )
    else:
        return 6 * ( r * ( 2 - 2 * a ) ) ** 2 * ( r + 4 * a ) / ( r_c ** 3 * ( r_c ** 3 - 18 * a ** 2 * r_c + 32 * a ** 3 ) )


def CollectDistances( structure, res_a, atom_a, res_b, atom_b, sequence_separation = 5, cutoff = 20 ):
    """
        @summary: Given a structure and a residue pair,
    """
    # collect residue pair distances
    a_b_distances = []
    for res_i in structure.get_residues():
        for res_j in structure.get_residues():
            if ( res_i.get_resname() == res_a and
                 res_j.get_resname() == res_b and
                 abs( res_i.id[1] - res_j.id[1] ) >= sequence_separation ):
                distance = res_i[atom_a] - res_j[atom_b]
                if distance <= cutoff:
                    a_b_distances.append( distance )
    return a_b_distances

def ParseCommandLineOptions():
    # parse command line options and arguments
    parser = ArgumentParser( description = "Collect the statistics of residue-residue contacts." )
    parser.add_argument( "-l", "--list", dest = "pdb_list", required = True, help = "PDB file" )
    parser.add_argument( "-p", "--atom_pair", dest = "atom_pair", required = True,
                             nargs = 4, help = "atom pair between which to measure the distance" )
    parser.add_argument( "-c", "--cutoff", dest = "cutoff", required = False, default = 20.0, type = float )
    parser.add_argument( "-s", "--sequence_separation", dest = "sequence_separation", required = False,
                        default = 6, type = int, help = "sequence separation beyond which residues will be considered" )
    parser.add_argument( "-o", "--outfile", dest = "outfile", required = True,
                         help = "output filename for contact statistics" )
    parser.add_argument( "-v", "--verbose", dest = "verbose", action = "store_true", help = "verbose mode" )
    return parser.parse_args()

if __name__ == "__main__":

    args = ParseCommandLineOptions()
    # print the collected arguments
    if args.verbose:
        print( "database:            " + args.pdb_list )
        print( "cutoff:              " + str( args.cutoff ) )
        print( "atom_pair:           " + str( args.atom_pair ) )
        print( "sequence_separation: " + str( args.sequence_separation ) )
        print( "outfile:             " + args.outfile )
        print( "verbose:             " + str( args.verbose ) )

    # collects statistics about residue pair distances
    a_b_distances = []
    res_a, atom_a, res_b, atom_b = args.atom_pair
    with open( args.pdb_list, "rt" ) as f:
        lines = f.readlines()
    for line in lines:
        pdb = line.rstrip()
        # parse the PDB file
        pdb_parser = PDBParser( PERMISSIVE = 1 )
        structure = pdb_parser.get_structure( str( pdb ), pdb )
        a_b_distances.extend( 
            CollectDistances( structure, res_a, atom_a, res_b, atom_b,
                              args.sequence_separation, args.cutoff
            )
        )

    bins = list( np.arange( 0, 15, 0.5 ) ) + [max( a_b_distances )]
    distance_pdf = np.histogram( a_b_distances, bins, normed = True )[0]
    normalization_factors = np.array( [ComputeNormalizationFunction( r, 15, 20 ) for r in bins[:-1]] )
    potentials = [10 if c else np.log( x ) for ( c, x ) in zip( distance_pdf == 0, distance_pdf / normalization_factors )]
    plt.plot( bins[1:], potentials )
    plt.savefig( args.outfile + 'pair_frequency.png' )
