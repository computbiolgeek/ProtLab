#!/usr/bin/env python3
"""
    @summary: This script takes a database in fasta format and a sequence id,
        it extracts the sequence corresponding to the given id from the database.
        the output is also in fasta format.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: 
    @copyright:     
"""

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# parse command line options and arguments
parser = ArgumentParser()
parser.add_argument("--db", help = "path to the sequence database in fasta format")
parser.add_argument("--id", help = "fasta id of the sequence to be extracted")
parser.add_argument("--out", help = "output filename")
args = parser.parse_args()

has_find = False
for record in SeqIO.parse(args.db, "fasta"):
    if record.id == args.id:
        print("Find sequence record for sequence id:", args.id)
        has_find = True
        SeqIO.write(record, args.out, "fasta")
        break

if not has_find:
    print("No sequence record for sequence id:", args.id, 
          "was found in database", args.db)
