"""
    @summary: 
"""

from argparse import ArgumentParser

# parse command line arguments
parser = ArgumentParser( description = """Computes position-specific conservation scores 
    given a multiple sequence alignment.""" )
parser.add_argument( '-i', '--input', dest = 'msa', required = True,
                    help = 'file containing the multiple sequence alignment' )
parser.add_argument( '-f', '--format', dest = 'format', required = True, default = 'fasta',
                    choices = ['fasta', 'clustal'], help = 'format of the multiple sequence alignment' )
parser.add_argument( '-t', '--type', dest = 'type', required = False, default = 'Protein',
                     choices = ['Protein', 'DNA', 'RNA'], help = 'type of sequences: DNA, RNA, Protein' )
parser.add_argument( '-o', '--output', dest = 'output', required = True,
                    help = 'output filename' )

args = parser.parse_args()


def read_fasta( filename ):
    """
        
    """
    # only sequences will be read in, ids are not of interest to this application
    seqs = []
    with open( filename, 'rt' ) as f:
        cur_seq = ''
        for line in f.readlines():
            if line.startswith( '>' ):
                # if cur_seq is not empty, it contains the preceding sequence
                if cur_seq != '':
                    seqs.append( cur_seq.upper() )
                    # set cur_seq to empty for reading the current sequence in the else block
                    cur_seq = ''
            else:
                cur_seq += line.strip()
        # add the last sequence
        seqs.append( cur_seq.upper() )
    return seqs
