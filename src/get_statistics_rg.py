"""
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu or comput.biol.geek@gmail.com
"""

from argparse import ArgumentParser
from pdb_utils import ReadPDB, ResHydrophobic, RadiusOfGyration
import csv
import sys

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

# compute radius of gyrations for each PDB
headers = ["filename", "number of residues", "rg_all", "rg_hydrophobic", "rg_ratio"]
records = []
with open( args.pdb_list, "rt" ) as f:
    for pdb in f.read().splitlines():
        residues, coordinates = ReadPDB( pdb )
        if not residues:
            print( pdb + ": no residues found, check the PDB file!" )
            sys.exit()
        # collect radius of gyrations
        is_hydrophobic = ResHydrophobic( residues )
        rg_all = RadiusOfGyration( coordinates )
        rg_hydrophobic = RadiusOfGyration( coordinates[is_hydrophobic] )
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
