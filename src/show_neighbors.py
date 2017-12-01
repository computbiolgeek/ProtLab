#!/usr/bin/env python3

import sys
from pymol import cmd

def main():
    '''
    '''
    args = sys.argv
    
    # if -h is given on the command line, print help
    if '-h' in args:
        print('usage: show_neighbors.py <pdb file> <residue ID>')
        sys.exit(0)

    # parse given PDB file and residue ID
    cmd.load(args[1])
    resi = int(args[2])
    
    # show cartoon
    cmd.show(representation='cartoon')
    cmd.cartoon('oval')
    cmd.set('cartoon_oval_length', value=1.0)
    # make selection
    seln = 'resi ' + resi
    cmd.select('res' + resi, selection=seln)
    cmd.show(representation='sphere', selection=seln)
    # select residues around resi 
    cmd.select()