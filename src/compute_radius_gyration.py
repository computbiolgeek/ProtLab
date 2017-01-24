"""
    @summary:
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu or comput.biol.geek@gmail.com
    @change:
    @copyright:
"""

from argparse import ArgumentParser
from collections import defaultdict
import csv
import sys

import numpy as np


# parse command line options and arguments
parser = ArgumentParser( description = "Compute the radius of gyration of the given \
protein and the radius of gyration of its hydrophobic residues." )
parser.add_argument( "-l", "--list", dest = "pdb_list", required = True, help = "PDB file" )
parser.add_argument( "-o", "--output", dest = "output", required = True, help = "output filename" )
parser.add_argument( "-v", "--verbose", dest = "verbose", help = "verbose mode" )
args = parser.parse_args()

# print the collected arguments
print( args.pdb_list )
print( args.output )
print( args.verbose )


def read_pdb(pdb_file):
    """
        @summary: Takes as input a PDB file, extracts residue names and the
            corresponding alpha carbon coordinates.
        @param pdb_file: the string filename of a PDB file
        @return:
            residue names: a length-N list containing the three-letter code for
                the type of each of the corresponding amino acid residues
            coordinates: a (N, 3) dimensional NumPy array containing the position
                of the alpha carbon in each of the N amino acid residues
    """
    res_names = []
    coordinates = []
    with open( pdb_file, "rt" ) as f:
        for line in f:
            if line.startswith( "ATOM" ):
                if line[12:16] == " CA ":
                    res_names.append( line[17:20] )
                    coordinates.append( [float( line[30:38] ),
                                              float( line[38:46] ),
                                              float( line[46:54] )] )
            elif line.startswith( "TER" ):
                break
    return res_names, np.array( coordinates )


def is_res_hydrophobic(residues):
    """
        @summary: Takes a length-N list of three-letter residue codes, returns
            a length-N array with Boolean values indicating whether the corresponding
            residue is hydrophobic or not.
        @param residues: a length-N list of three-letter residue codes
        @return:
            is_hydrophobic, a length-N NumPy array with True or False boolean
            values indicating whether the residue codes in "residues" correspond
            hydrophobic or hydrophilic residues
    """
    hydrophobic = {"ALA",
                   "CYS",
                   "PHE",
                   "ILE",
                   "LEU",
                   "MET",
                   "VAL"}
    is_hydrophobic = [residue in hydrophobic for residue in residues]
    return np.array( is_hydrophobic )


def radius_of_gyration(coordinates):
    """
        @summary: Computes the radius of gyration of the passed set of atoms.
        @param
            coordinates: a (N, 3) dimensional NumPy array.
        @return:
            Rg: a float with the corresponding radius of gyration
    """
    r_mean = coordinates.mean( axis = 0 )
    rg = np.mean( np.sum( ( coordinates - r_mean ) ** 2 ) )
    return np.sqrt( rg )


def get_contacts(coordinates, residues, cutoff):
    """
        @summary:
        @param
            coordinates:
            residues:
            cutoff:
        @return:
            res_freqs:
            contacts:
    """
    amino_acids = ["ALA",
                   "CYS",
                   "ASP",
                   "GLU",
                   "PHE",
                   "GLY",
                   "HIS",
                   "ILE",
                   "LYS",
                   "LEU",
                   "MET",
                   "ASN",
                   "PRO",
                   "GLN",
                   "ARG",
                   "SER",
                   "THR",
                   "VAL",
                   "TRP",
                   "TYR"]
    res_freqs = {aa: residues.count( aa ) for aa in amino_acids}
    contacts = defaultdict( list )
    # number of residues
    n = coordinates.shape[0]
    # compute the upper triangular matrix
    for i in range( n ):
        for j in range( i + 1, n ):
            dist_ij = np.sqrt( ( sum( coordinates[i] - coordinates[j] ) ** 2 ) )
            if dist_ij < cutoff:
                contacts[( residues[i], residues[j] )].append( dist_ij )
    return res_freqs, contacts

# compute radius of gyrations for each PDB
headers = ["filename", "number of residues", "rg_all", "rg_hydrophobic", "rg_ratio"]
records = []
contacts_all = defaultdict( list )
with open( args.pdb_list, "rt" ) as f:
    for pdb in f.read().splitlines():
        residues, coordinates = read_pdb(pdb)
        if not residues:
            print( pdb + ": no residues found, check the PDB file!" )
            sys.exit()
        # collect all contacts
        contacts = get_contacts(coordinates, residues, 9)
        for key, value in contacts.items():
            contacts_all[key].extend( value )
        # collect radius of gyrations
        is_hydrophobic = is_res_hydrophobic(residues)
        rg_all = radius_of_gyration(coordinates)
        rg_hydrophobic = radius_of_gyration(coordinates[is_hydrophobic])
        rg_ratio = rg_hydrophobic / rg_all
        record_current = {headers[0]: pdb,
                          headers[1]: len( residues ),
                          headers[2]: rg_all,
                          headers[3]: rg_hydrophobic,
                          headers[4]: rg_ratio}
        records.append( record_current )
# print a summary on the radius of gyration of the given protein
with open( args.output, "wt" ) as csv_file:
    csv_writer = csv.DictWriter( csv_file, headers )
    csv_writer.writeheader()
    csv_writer.writerows( records )

print( contacts_all[( "GLU", "LYS" )] )
