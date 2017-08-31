#!/usr/bin/env python3

from Bio import pairwise2, Seq, SeqIO, SeqRecord, Align, AlignIO
from Bio.SubsMat import MatrixInfo
from Bio.Alphabet import IUPAC
from argparse import ArgumentParser
import sys

def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input')
    parser.add_argument('-o', '--output', dest='output')
    args = parser.parse_args()
    
    records = list(SeqIO.parse(args.input, 'fasta', alphabet=IUPAC.protein))
    if len(records) != 2:
        print('Fasta file has ' + str(len(records)) + ' sequences. It must contain two sequences')
        sys.exit(0)
    # print the given two sequences
    print(records[0], '\n')
    print(records[1], '\n')
    
    # now align the two given sequences
    print('Now aligning the two given sequences using global alignment and BLOSUM62 substitution matrix ...')
    alignments = pairwise2.align.globaldx(records[0].seq, records[1].seq, MatrixInfo.blosum62)
    multiple_alignment = Align.MultipleSeqAlignment([
            SeqRecord.SeqRecord(seq=Seq.Seq(alignments[0][0]), id=records[0].id, description=records[0].description),
            SeqRecord.SeqRecord(seq=Seq.Seq(alignments[0][1]), id=records[1].id, description=records[1].description)
            ])
    AlignIO.write(multiple_alignment, args.output, 'fasta')
    print('The following alignment has been written to ' + args.output, '\n')
    print(pairwise2.format_alignment(*alignments[0]))

if __name__ == '__main__':
    main()