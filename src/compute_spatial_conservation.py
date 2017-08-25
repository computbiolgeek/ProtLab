"""
    @summary: This program computes residue spatial conservation index which is defined
        as the sum of c/r^2 across all neighbor residues. c is the conservation of a
        neighboring residue and r is the distance between the residue of interest and the
        neighbor residue. This index roughly measures how spatially constrained in terms of
        evolutionary conservation.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @copyright: Bian Li
    @change: 12/14/16 implemented the first version of the algorithm (version 0.10)
"""

from argparse import ArgumentParser
from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import csv

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

class NeighborList( object ):
    """
        @summary: A wrapper for a list of Biopython Residue objects.
    """
    def __init__( self, neighbor_list ):
        self._neighbor_list = neighbor_list

    def Count( self, residue ):
        """
            @summary: Counts the number of 'residue' in the neighbor list.
        """
        count = 0
        for res in self._neighbor_list:
            if res.get_resname() == residue:
                count += 1
        return count;

    def Length( self ):
        return len( self._neighbor_list )

    def Append( self, new_neighbor ):
        self._neighbor_list.append( new_neighbor )

    def GetAllCoords( self, atom ):
        """
            @summary: Returns a numpy array of coordinates of all residues in the neighbor list.
            @param atom: one of {"CA", "CB", "Centroid"}
        """
        if atom == "Centroid":
            return np.array( [ComputeSideChainCentroid( res ) for res in self._neighbor_list] )
        else:
            return np.array( [ res["CA"].get_coord() if res.get_resname() == "GLY"
                              else res[atom].get_coord() for res in self._neighbor_list] )

    def GetAllBFactors( self ):
        """
            @summary: Returns a numpy array of b factors of all residues in the neighbor list.
        """
        return np.array( [res["CA"].get_bfactor() for res in self._neighbor_list] )


def ParseBFactors( b_factor_file ):
    """
        @summary: Parses a .csv b-factor file into a dictionary.
    """
    with open( b_factor_file, "rt" ) as f:
        reader = csv.reader( f )
        d = {int( row[0] ): float( row[1] ) for row in reader if row}
    return d

def SetBFactors( chain, b_factors ):
    """
        @summary: Sets residue b-factors to given values.

        @param chain: a Biopython Chain object
        @param b_factors: a dictionary which has residue IDs as keys
            and given b-factor values as values
    """
    for res in chain:
        res_id = res.get_id()[1]
        for atom in res.get_list():
            atom.set_bfactor( b_factors[res_id] )


def ComputeSpatialConservation( residue, atom , neighbor_list ):
    """
        @summary: Computes residue spatial conservation.

        @param residue: a Biopython Residue object for which to compute spatial conservation
        @param atom: atom type used to measure distance
        @param neighbor_list: an object of type NeighborList
    """
    assert neighbor_list.Length() != 0
    neighbor_p = neighbor_list.GetAllCoords( atom )
    if atom == "Centroid":
        residue_p = ComputeSideChainCentroid( residue )
    elif residue.get_resname() == "GLY":
        residue_p = residue["CA"].get_coord()
    else:
        residue_p = residue[atom].get_coord()
    d_ijs = np.linalg.norm( residue_p - neighbor_p, axis = 1 )
    # get all b factors
    b_factors = neighbor_list.GetAllBFactors()
    return np.sum( b_factors * ( 1.0 / d_ijs ** 2 ) )

def main():
    """
    """
    parser = ArgumentParser( description = "Hi there, I compute residue spatial conservation for you!" )
    parser.add_argument( "-p", "--pdb", dest = "pdb", required = True,
                        help = "PDB file of mono-chain protein or homo-oligomers" )
    parser.add_argument( "-b", "--b_factors", dest = "b_factors", required = True,
                         help = "file that stores residue conservations computed from multiple sequence alignment" )
    parser.add_argument( "-o", "--outfile", dest = "outfile", required = True,
                        help = "output file name" )
    parser.add_argument( "-a", "--atom", dest = "atom", choices = {"CA", "CB", "Centroid"},
                        default = "Centroid", help = "atom type used to measure distance" )
    parser.add_argument( "-r", "--radius", dest = "radius", type = float, default = 10.0,
                         help = "radius for computing neighor list" )
    parser.add_argument( "-s", "--sequence_separation", dest = "sequence_separation", type = int,
                         default = 3, help = "sequence separation beyond which "
                        "residues will be considered" )
    parser.add_argument( "-v", "--verbose", dest = "verbose", action = "store_true",
                        help = "verbose mode" )
    args = parser.parse_args()

    # print collected command line arguments
    if args.verbose:
        print( "pdb:       " + args.pdb )
        print( "b_factors: " + args.b_factors )
        print( "outfile:   " + args.outfile )
        print( "atom:      " + args.atom )
        print( "radius:    " + str( args.radius ) )
        print( "verbose:   " + str( args.verbose ) )

    # parse the given PDB file
    pdb_parser = PDBParser( PERMISSIVE = 1 )
    structure = pdb_parser.get_structure( "given_pdb", args.pdb )

    # parse the given b factors into a dictionary
    bfactor_dict = ParseBFactors( args.b_factors )

    # set b factors (residue conservations) for all residues
    for model in structure.get_list():
        for chain in model.get_list():
            SetBFactors( chain, bfactor_dict )

    # compute contact number for each residue
    atom_list = list( structure.get_atoms() )
    ns = NeighborSearch( atom_list )
    spatial_conservations = []
    for model in structure.get_list():
        for chain in model.get_list():
            # chain_id = chain.get_id()
            for residue in chain.get_list():
                if args.atom == "Centroid":
                    measurement_point = ComputeSideChainCentroid( residue )
                elif residue.get_resname() == "GLY":
                    measurement_point = residue["CA"].get_coord()
                else:
                    measurement_point = residue[args.atom].get_coord()
                # compile a list of neighbors for the current residue
                neighbor_list = [res for res in ns.search( measurement_point, args.radius, 'R' )
                                 if res.get_id()[1] - residue.get_id()[1] > args.sequence_separation
                                    or res.get_id()[1] - residue.get_id()[1] < -args.sequence_separation]
                # if neighbor list is empty, set spatial conservation to 0
                if not neighbor_list:
                    current_res_sc = 0
                else:
                    current_res_sc = ComputeSpatialConservation( residue, args.atom, NeighborList( neighbor_list ) )
                residue_id = residue.get_id()[1]
                spatial_conservations.append( 
                    ( 
                      residue_id,
                      current_res_sc
                    )
                )

    # write contact numbers into a .csv file
    with open( args.outfile, "wt" ) as f:
        csv_f = csv.writer( f )
        for row in spatial_conservations:
            csv_f.writerow( row )

if __name__ == "__main__":
    main()
