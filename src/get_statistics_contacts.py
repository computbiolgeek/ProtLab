"""
    @summary: 
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu or comput.biol.geek@gmail.com
    @change: 
    @copyright: 
"""

from argparse import ArgumentParser
import numpy as np
from collections import defaultdict
from pdb_utils import *
import csv
import sys
import itertools

# parse command line options and arguments
parser = ArgumentParser(description="Collect the statistics of residue-residue contacts.")
parser.add_argument("-l", "--list", dest="pdb_list", required=True, help="PDB file")
parser.add_argument("-c", "--cutoff", dest="cutoff", required=False, default=8.0, type=float)
parser.add_argument("-f", "--aa_frequency", dest="frequency", required=False, 
                    default="aa_frequency.csv", help="output filename for residue frequencies")
parser.add_argument("-s", "--contact_statistics", dest="contacts", required=False,
                    default="contact_statistics.csv", help="output filename for contact statistics" )
parser.add_argument("-v", "--verbose", dest="verbose", help="verbose mode")
args = parser.parse_args()

# print the collected arguments
print(args.pdb_list)

# collects statistics about contacts
contacts_all = defaultdict(list)
frequencies_all = defaultdict(int)
with open(args.pdb_list, "rt") as f:
    for pdb in f.read().splitlines():
        residues, coordinates = ReadPDB(pdb, atom="CB")
        if not residues:
            print(pdb + ": no residues found, check the PDB file!")
            sys.exit()
        # collect all contacts
        residue_frequencies, contacts = GetContacts(coordinates, residues, args.cutoff)
        for key, value in residue_frequencies.items():
            frequencies_all[key] += value
        for key, value in contacts.items():
            contacts_all[key].extend(value)

# output statistics about contacts
with open(args.contacts, "wt") as f:
    csv_writer = csv.writer(f)
    csv_writer.writerow(contacts_all.keys())
    csv_writer.writerows(itertools.izip_longest(*contacts_all.values()))