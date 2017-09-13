#!/usr/bin/env python

# import required modules
from argparse import ArgumentParser
import numpy as np
import csv

# create a parser and parse command line arguments
parser = ArgumentParser()
parser.add_argument('-f', '--fasta', dest='fasta', help='fasta file which contains the sequence of residues '
                    'that the OPM pdb file has resolved coordiates.')
parser.add_argument('-r', '--resseq', dest='resseq', help='file that contains residue resseq ids, extracted from OPM pdb file')
parser.add_argument('-s', '--segments', dest='segments', help='segment specifitions provided by OPM')
parser.add_argument('-o', '--output', dest='output', help='output file that contains the mapping between residues and TM segment annotations')
args = parser.parse_args()

with open(args.segments, 'rt') as f:
    for l in f.readlines():
        fields = l.strip().split(',')
        pdb_id = fields[0]
        segments = []
        for k in fields[1:]:
            # add transmembrane segments sequence ranges
            segments.append(tuple(k.split(':')[1].strip().split('-')))
            # cast strings to integers
            segments = [(int(r[0]), int(r[1])) for r in segments]

# read in sequence
with open(args.fasta, 'rt') as f:
    lines = [l.strip() for l in f.readlines() if not l.startswith('>')]
    sequence = ''.join(lines)

# annotation
annotation = np.array(['-' for i in range(len(sequence))])
for seg in segments:
    annotation[(seg[0]-1):seg[1]] = 'M'
annotation = ''.join(annotation)

# write annotation along with amino acid positions and names into a csv file
with open(args.output, 'wt') as f:
    csv_writer = csv.writer(f)
    i = 1
    for s, a in zip(sequence, annotation):
        csv_writer.writerow((i, s, a))
        i += 1
