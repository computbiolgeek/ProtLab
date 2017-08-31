#!/usr/bin/env python3

# import required modules
from Bio import PDB, SeqIO, SeqRecord, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from argparse import ArgumentParser
from os.path import basename


def main():
    # parse command line arguments
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', help='input PDB file')
    parser.add_argument('-c', '--chains', dest='chains', help='chains to extract sequence for')
    parser.add_argument('-o', '--output', dest='output', help='output fasta file')
    args = parser.parse_args()
    
    print( "input file:          " + args.input )
    print( "output file:         " + args.output )
    print( "chains:              " + args.chains )
    
    # extract sequences from SEQRES records
    pdb_id = basename(args.input).split('.')[0]
    with open(args.input, 'rt') as f:
        seqres_sequences = list(SeqIO.parse(f, 'pdb-seqres'))
    
    # extract sequences from residues with resolved coordinates
    structure = PDB.PDBParser().get_structure(id=pdb_id, file=args.input)
    model = structure[0]
    peptide_builder = PDB.Polypeptide.PPBuilder()
    sequences = []
    for c in args.chains:
        chain = peptide_builder.build_peptides(model[c])
        coord_sequence = chain[0].get_sequence()
        print('Sequence for chain ' + c + ' extracted from coordinates:')
        print(coord_sequence, '\n')
        # if the SEQRES records are missing from the PDB file, use sequence from coordinates
        if len(seqres_sequences) == 0:
            print('SEQRES records in the given PDB file are missing! Using sequence'
                  'extracted from coordinates.')
            sequence = coord_sequence
        else:
            for record in seqres_sequences:
                if c == record.id:
                    print('Sequence for chain ' + c + ' in SEQRES:')
                    print(record.seq, '\n')
            # pairwise alignment between the two sequences
            print('Here is an alignment of the two sequences:')
            alignment = pairwise2.align.globaldx(coord_sequence, record.seq, matlist.blosum62)
            print(pairwise2.format_alignment(*alignment[0]))
            
            # ask for which sequence to choose
            s = input('Which sequence are you interested? 1 for sequence from coordinates '
                      '2 for sequence from SEQRES: ')
            if int(s) == 1:
                sequence = coord_sequence
            else:
                sequence = record.seq
            
        # append sequence for the current chain
        sequences.append(
                SeqRecord.SeqRecord(
                        seq=sequence, 
                        id='chain ' + c, 
                        description='sequence for chain '+ c + ' from ' + pdb_id))
    
    # write sequences to a fasta file
    with open(args.output, 'wt') as f:
        SeqIO.write(sequences, f, 'fasta')
        
if __name__ == '__main__':
    main()