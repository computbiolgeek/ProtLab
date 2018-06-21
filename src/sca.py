#!/usr/bin/env python3

import sys
from Bio.AlignIO import AlignIO


def create_msa_array(msa):
    '''
        remove columns where the subject sequence has gaps, the top sequence is assumed
        to be the subject sequence and gaps are assumed to be represented as '-'s
    '''
    print('Computing a NumPy ndarray representation of the alignment')
    msa_arr = np.asarray([list(seq.upper()) for seq in msa])
    msa_arr[msa_arr == 'B'] = 'D'
    msa_arr[msa_arr == 'Z'] = 'Q'
    return msa_arr


def main():
    '''
    Statistical coupling analysis
    
