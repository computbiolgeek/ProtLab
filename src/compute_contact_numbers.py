#!/usr/bin/env python

"""
    @summary:
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: 12/12/16 the argument to the np.cos() function was converted to radians
             25/08/17 changed get_atom() function call to get_list() for compatibility with older Biopython versions.
                      added an option for considering residues only from different chains. 
             
"""

from argparse import ArgumentParser
import csv
from Bio.PDB import PDBParser, NeighborSearch
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
    for atom in residue.get_list():
        if atom.get_name() not in backbone_atoms:
            centroid += atom.get_coord()
            number_sidechain_atoms += 1
    return centroid / number_sidechain_atoms

def ComputeContactNumber( residue,
                          neighbor_list,
                          bounds,
                          measurement_point = "Centroid",
                          cn_type = "unweighted"
                          ):
    """
        @summary: Computes the contact number for the given residue and its neighbor list.

        @param residue: amino acid residue for which to compute the contact number
        @param neighbor_list: a list of neighboring amino acid resdiues
        @param bounds: a list containing the lower distance and upper distance bounds
        @param measurement_point: one of {"CA", "CB", "Centroid"}
        @param cn_type: weighted or unweighted

        @see: Bian Li, et al JCIM 2016 for algorithmic details

    """
    # compute contact number
    cn = 0
    d_ijs = []
    if measurement_point == "Centroid":
        residue_p = ComputeSideChainCentroid( residue )
        for neighbor in neighbor_list:
            neighbor_p = ComputeSideChainCentroid( neighbor )
            d_ijs.append( np.linalg.norm( residue_p - neighbor_p ) )
    else:
        if residue.get_resname() == "GLY":
            residue_p = residue["CA"].get_coord()
        else:
            residue_p = residue[measurement_point].get_coord()
        for neighbor in neighbor_list:
            if neighbor.get_resname() == "GLY":
                neighbor_p = neighbor["CA"].get_coord()
            else:
                neighbor_p = neighbor[measurement_point].get_coord()
            d_ijs.append( np.linalg.norm( residue_p - neighbor_p ) )
    for d_ij in d_ijs:
        if d_ij <= bounds[0]:
            cn += 1
        elif d_ij > bounds[0] and d_ij <= bounds[1]:
            if cn_type == "unweighted":
                cn += 1
            else:
                cn += ( 1.0 / 2 ) * ( 
                    np.cos( ( d_ij - bounds[0] ) / ( bounds[1] - bounds[0] ) * np.pi ) + 1
                )
    return cn

def main():
    """
    """
    # specify command line flags
    parser = ArgumentParser( description = "Hi there, I compute contact numbers for you." )
    parser.add_argument( "-i", "--infile", dest = "infile", required = True,
                        help = "input PDB file for computing contact numbers" )
    parser.add_argument( "-o", "--outfile", dest = "outfile", required = True,
                        help = "output file to store computed contact numbers" )
    parser.add_argument( "-t", "--type", dest = "type", choices = {"weighted", "unweighted"},
                         default = "weighted", help = "the type of contact number definition" )
    parser.add_argument( "-b", "--bounds", dest = "bounds", nargs = 2, type = float,
                         default = [4, 11.4], help = "cutoff distances within which residues will be considered" )
    parser.add_argument( "-s", "--sequence_separation", dest = "sequence_separation", type = int,
                         default = 3, help = "sequence separation beyond which "
                        "residues will be considered" )
    parser.add_argument( "-c", "--different_chains", dest = "different_chains", action = "store_true", default = False, 
                        help = "true if consider only residues from different chains, else consider all chains")
    parser.add_argument( "-m", "--measurement_point", dest = "measurement_point",
                         choices = {"CA", "CB", "Centroid"}, default = "Centroid",
                         help = "points between which to measure the distance" )
    parser.add_argument( "-v", "--verbose", dest = "verbose", action = "store_true",
                        help = "verbose mode" )
    # parse command line arguments
    args = parser.parse_args()

    # print collected command line arguments
    if args.verbose:
        print( "input file:          " + args.infile )
        print( "output file:         " + args.outfile )
        print( "type:                " + args.type )
        print( "bounds:              " + str( args.bounds ) )
        print( "sequence_separation: " + str( args.sequence_separation ) )
        print( "measurement_point:   " + args.measurement_point )
        print( "different chains:    " + str( args.different_chains ) )
        print( "verbose:             " + str( args.verbose ) )

    # parse the given PDB file
    pdb_parser = PDBParser( PERMISSIVE = 1 )
    structure = pdb_parser.get_structure( "pdb", args.infile )

    # compute contact number for each residue
    model = structure[0] # consider only the first model if the PDB file has multiple models 
    atom_list = list( model.get_atoms() )
    ns = NeighborSearch( atom_list )
    contact_numbers = []
    for residue in model.get_residues():
        if args.measurement_point == "CA" or residue.get_resname() == "GLY":
            measurement_point = residue["CA"].get_coord()
        elif args.measurement_point == "CB":
                measurement_point = residue["CB"].get_coord()
        else:
            measurement_point = ComputeSideChainCentroid( residue )
        # compile a list of neighbors for the current residue
        neighbor_list = list(ns.search( measurement_point, args.bounds[1], 'R' ))
        # if consider only different chains
        if args.different_chains:
            neighbor_list = [res for res in neighbor_list if res.get_full_id()[2] != residue.get_full_id()[2]]
        else:
            neighbor_list = [res for res in neighbor_list
                         if res.get_id()[1] - residue.get_id()[1] > args.sequence_separation
                            or res.get_id()[1] - residue.get_id()[1] < -args.sequence_separation]
        residue_id = residue.get_id()[1]
        chain_id = residue.get_full_id()[2]
        contact_numbers.append( 
            ( residue_id, chain_id, ComputeContactNumber( residue, neighbor_list, args.bounds,
                                                args.measurement_point, args.type ) )
        )

    # write contact numbers into a .csv file
    with open( args.outfile, "wt" ) as f:
        csv_f = csv.writer( f )
        for row in contact_numbers:
            csv_f.writerow( row )

if __name__ == "__main__":
    main()
