#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
from Bio import AlignIO
import re
import sys
import math
import pandas as pd


# alphabet, global variable
alphabet = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-']


def create_msa_array(msa):
    '''
        remove columns where the subject sequence has gaps, the top sequence is assumed
        to be the subject sequence and gaps are assumed to be represented as '-'s
    '''
    print('Computing a NumPy ndarray representation of the alignment')
    msa_arr = np.array([list(rec.upper()) for rec in msa])
    msa_arr[msa_arr == 'B'] = 'D'
    msa_arr[msa_arr == 'Z'] = 'Q'
    return msa_arr


def compute_aa_probs(column_i):
    '''
    '''
    # amino acid background relative frequencies 
    # bg_rf = np.array([0.085786, 0.045676, 0.047306, 0.058022, 0.018036, 0.037722, 0.059724, 0.081155, 0.021639, 
    #        0.052944, 0.081156, 0.058717, 0.021109, 0.039946, 0.048178,0.063047, 0.060835, 0.014256, 0.036310, 0.068436])

    f = np.array([float(sum(column_i == aa) + 1) for aa in alphabet])
    # print('The frequency of each amino acid type at position ' + str(i) + ' is: ')
    # print(f)

    # compute relative frequencies
    rf = np.true_divide(f, sum(f))
    # print('The probability that each amino acid type at position ' + str(i) + ' is: ')
    # print(rf)

    return rf


def compute_joint_probs(column_i, column_j):
    '''
    Compute the joint distribution of amino acid pairs for position pair i, j.
    See Weigt et al, PNAS, 2009 for technical details.
    '''
    if len(column_i) != len(column_j):
        print('Unequal depth of the two given alignment column pair.')
        sys.exit(0)
    joint_probs = []
    for a in alphabet:
        # count frequency of amino acid pair (a, b) in columns (i, j)
        joint_freq = np.array([sum((column_i == a) * (column_j == b)) + (1. / 20) for b in alphabet])

        # compute joint probability
        joint_prob = joint_freq / (20 + len(column_i))
        joint_probs.append(joint_prob)

    return np.array(joint_probs)


def compute_mi(probs_i, probs_j, joint_probs_ij):
    '''
    Compute the mutual information between position pair(i, j) from the given multiple sequence alignment.
    See Weigt et al, PNAS, 2009 for technical details.
    '''
    mi = 0
    for a in range(21):
        for b in range(21):
            mi += joint_probs_ij[a, b] * math.log(joint_probs_ij[a, b] / (probs_i[a] * probs_j[b]))
    return mi


def compute_excess_mi(column_i, column_j):
    '''
    '''
    probs_i = compute_aa_probs(column_i)
    probs_j = compute_aa_probs(column_j)
    joint_probs = compute_joint_probs(column_i, column_j)
    apparent_mi = compute_mi(probs_i, probs_j, joint_probs)
    print('The apparent mutual information is:', apparent_mi)
    n_permutations = 1000
    empirical_mis = []
    for i in range(n_permutations):
        column_j_permuted = np.random.permutation(column_j)
        joint_probs = compute_joint_probs(column_i, column_j_permuted)
        mi_permutation = compute_mi(probs_i, probs_j, joint_probs)
        print('Mutual information from permutation round', i, 'is:', mi_permutation)
        empirical_mis.append(mi_permutation)
    n = sum(np.array(empirical_mis) > apparent_mi)
    if n >= 10:
        print('Got', n, 'permutations where the mutual information is greater than the apparent mutual information:', apparent_mi)
        return np.NAN
    else:
        return apparent_mi - np.mean(empirical_mis)


def main():
    '''
    '''
    # parser for parsing command line arguments
    parser = ArgumentParser(description='Compute mutual information from a multiple sequence alignment for a given list of pairs of residues.')
    parser.add_argument('-m', '--msa', dest='msa', help='multiple sequence alignment in fasta format')
    parser.add_argument('-r', '--res', dest='res', nargs=2, help='two files, each containing a list of residue numbers')
    parser.add_argument('-o', '--output', dest='output', help='a matrix of pairwise mutual information')

    # parse command line arguments
    args = parser.parse_args()

    # read the multiple sequence alignment
    alignment = AlignIO.read(args.msa, format='fasta')

    # create a 2D numpy array from the alignment
    msa_array = create_msa_array(alignment)
    # print(msa_array.shape)

    # residues given on the command line
    with open(args.res[0], 'rt') as f:
        res_ids_a = [tuple(re.split(r'[,;\s]\s*', line.strip())) for line in f.readlines()]
    residues_a = [int(i[1]) for i in res_ids_a]
    with open(args.res[1], 'rt') as f:
        res_ids_b = [tuple(re.split(r'[,;\s]\s*', line.strip())) for line in f.readlines()]
    residues_b = [int(i[1]) for i in res_ids_b]

    # compute distributions of amino acids at each position
    aa_probs_a = [compute_aa_probs(msa_array[:, i-1]) for i in residues_a]
    aa_probs_b = [compute_aa_probs(msa_array[:, i-1]) for i in residues_b]
    
    # compute the mutual information between the distributions of amino acids occurring in each pair of positions
    mi = []
    for i, a in enumerate(residues_a):
        mi_a = []
        if aa_probs_a[i][-1] > 0.5:
            mi_a = [np.NAN for i in range(len(residues_b))]
            mi.append(mi_a)
            continue
        for j, b in enumerate(residues_b):
            if aa_probs_b[j][-1] > 0.5:
                mi_a.append(np.NAN)
                continue
            print('Now computing the joint distribution of amino acids at position pair: ' + '(' + str(a) + ', ' + str(b) + ')')
            joint_probs_ab = compute_joint_probs(msa_array[:, a-1], msa_array[:, b-1])
            mi_a.append(compute_excess_mi(msa_array[:, a-1], msa_array[:, b-1]))
        mi.append(mi_a)

    # make a DataFrame for formatting output
    subject_seq = msa_array[0]
    indices = [''.join([i[0], i[1], subject_seq[int(i[1])-1]]) for i in res_ids_a]
    headers = [''.join([i[0], i[1], subject_seq[int(i[1])-1]]) for i in res_ids_b]
    mi_df = pd.DataFrame(np.array(mi), index=indices, columns=headers)

    # write to file
    mi_df.to_csv(args.output, index=True, header=True, sep=',', float_format='%.4f')


if __name__ == '__main__':
    main()
