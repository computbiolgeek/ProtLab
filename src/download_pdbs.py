"""
    @summary: This script takes in a list of PDB IDs and retrieves the PDB coordinate file for each ID.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: 1/17/2017
"""

from argparse import ArgumentParser
import os
from Bio.PDB import PDBList, PDBParser, PDBIO


# parse command line options and arguments
parser = ArgumentParser( description = "Hi there, I download PDB file for you in batch!" )
parser.add_argument( "-l", "--list", dest = "list",
                     help = "a file containing a list of PDB IDs one per line" )
parser.add_argument( "-c", "--chains", dest = "chains", action = "store_true",
                     help = "whether the PDB IDs contains chain IDs" )
parser.add_argument( "-d", "--dir", dest = "dir",
                     help = "path to the directory where to store the PDB files" )
parser.add_argument( "-v", "--verbose", dest = "verbose", action = "store_true",
                     help = "print details while downloading" )
args = parser.parse_args()

# print fetched arguments
if args.verbose:
    print( "list: " + args.list )
    print( "dir:  " + args.dir )

# read the PDB IDs into a python list
with open( args.list, "rt" ) as f:
    pdb_ids = f.read().splitlines()

# retrieve PDB files and extract essential headers
pdbl = PDBList()
for pdb_id in pdb_ids:
    four_letter = pdb_id[:4].lower()
    # retrieve current PDB file and store it under a directory tree
    if not os.path.exists( args.dir + "pdb" + four_letter + ".ent" ):
        try:
            pdbl.retrieve_pdb_file( four_letter, pdir = args.dir )
        except Exception:
            print( "Exception occurred while retrieving the PDB file for " + four_letter )
            continue
    # if given PDB IDs are actually chains, then get the chain
    if args.chains:
        chain_id = pdb_id[-1]
        chain_file = args.dir + four_letter + chain_id + ".pdb"
        if os.path.exists( chain_file ):
            print( chain_file + " already exists, skip ..." )
        else:
            pdb_parser = PDBParser( PERMISSIVE = 1 )
            structure = pdb_parser.get_structure( four_letter, args.dir + "pdb" + four_letter + ".ent" )
            chain = structure[0][chain_id]
            pdb_writer = PDBIO()
            pdb_writer.set_structure( chain )
            pdb_writer.save( chain_file )
