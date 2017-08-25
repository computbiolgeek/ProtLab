#!/usr/bin/python
"""
    @summary: The script gives three essential attributes (structure method, resolution, 
        and journal reference) of each given PDB ID.
    @author: Bian Li
    @change: Aug 02, 2016
    @copyright: Bian Li
    @contact: bian.li@vanderbilt.edu or libian1128@gmail.com
"""

from Bio.PDB import *
from argparse import ArgumentParser
import csv
import os

# parse command line options and arguments
parser = ArgumentParser()
parser.add_argument("--pdbs", help = "a file containing a list of PDB IDs one per line")
parser.add_argument("--output", help = "output file name, the format is \"comma separated values\"")
args = parser.parse_args()

# read the PDB IDs into a python list
with open(args.pdbs, "rt") as f:
    pdbs = f.read().splitlines()

# retrieve PDB files and extract essential headers 
headers = ["PDB ID", "Structure Method", "Resolution", "Journal Reference"]
all_records = []
pdbl = PDBList()
pdb_parser = PDBParser(PERMISSIVE=1)
for pdb in pdbs:
    pdbl.retrieve_pdb_file(pdb, pdir = "/tmp/")
    pdb_filename = "/tmp/pdb" + pdb + ".ent"
    structure = pdb_parser.get_structure(pdb, pdb_filename)
    current_record = {headers[0]: pdb, headers[1]: structure.header["structure_method"],
                      headers[2]: structure.header["resolution"], headers[3]: structure.header["journal_reference"]}
    all_records.append(current_record)
    os.remove(pdb_filename)

# write the results into the specified file
with open(args.output, "wt") as csvfile:
    csv_writer = csv.DictWriter(csvfile, headers)
    csv_writer.writeheader()
    csv_writer.writerows(all_records)