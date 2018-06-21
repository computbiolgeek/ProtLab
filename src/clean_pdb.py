#!/usr/bin/env python3

"""
    @summary: This is a simple script that does PDB file cleaning. By default, it does not split chains
        and everything other than amino acid residues are removed. However, chain splitting can be specified
        by using the flag --split_chains. In this case, protein chains are cleaned whereas non-protein chains
        are saved directly without cleaning.
    @author: Bian Li
    @copyright: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: 12/12/16
"""

from argparse import ArgumentParser
from Bio.PDB import PDBParser, Select, is_aa, PDBIO

parser = ArgumentParser( description = "I'll clean the PDB file you gave me!" )
parser.add_argument( "-i", "--infile", dest = "infile", required = True,
                     help = "PDB file to be cleaned" )
parser.add_argument( "-o", "--outfile_prefix", dest = "outfile_prefix", required = True,
                     help = "output filename" )
parser.add_argument( "-s", "--split_chains", dest = "split_chains", required = False,
                    action = "store_true", help = "whether to split chains or not for a multi-chain protein" )
parser.add_argument( "-v", "--verbose", dest = "verbose", required = False,
                     action = "store_true", help = "verbose mode" )
args = parser.parse_args()

# print the collected command line arguments
if args.verbose:
    print( args.infile )
    print( args.outfile_prefix )
    print( args.split_chains )
    print( args.verbose )


class AminoAcidSelect( Select ):
    """
        Subclassing Select while overriding accept_residue() method.
    """
    def accept_residue( self, residue ):
        if is_aa( residue ):
            return True
        else:
            return False

if __name__ == "__main__":
    # parse input PDB file, PERMISSIVE = 1 ignores a number of common errors in PDB files
    pdb_parser = PDBParser( PERMISSIVE = 1 )
    structure = pdb_parser.get_structure( "pdb", args.infile )
    # do PDB file cleaning
    pdb_io = PDBIO()
    if not args.split_chains:
        pdb_io.set_structure( structure )
        pdb_io.save( args.outfile_prefix + "_cleaned.pdb", AminoAcidSelect() )
    else:
        print( "Splitting chains ..." )
        for chain in structure.get_chains():
            pdb_io.set_structure( chain )
            chain_id = chain.get_id()
            print( "Found chain " + chain_id + " and let's see whether it is a protein chain or not" )
            has_aa = has_non_aa = False
            for residue in chain:
                if args.verbose:
                    print( "Processing " + str( residue ) )
                if is_aa( residue ):
                    has_aa = True
                else:
                    if args.verbose:
                        print( "Found " + str( residue ) + ", which is not an amino acid residue" )
                    has_non_aa = True
                if has_aa and has_non_aa:
                    break
            if not has_aa and has_non_aa:
                print( "This chain only contains stuff other than amino acid residues, it's not a protein chain! "
                "Save without cleaning!" )
                pdb_io.save( args.outfile_prefix + "_" + chain_id + ".pdb" )
            else:
                print( "This chain contains amino acid residues, presumably it's a protein chain, lets clean it!" )
                pdb_io.save( args.outfile_prefix + "_" + chain_id + "_cleaned.pdb", AminoAcidSelect() )

