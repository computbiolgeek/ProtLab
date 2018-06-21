#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
from Bio import AlignIO
from Bio.PDB import PDBParser
import re
import sys
import math
import pandas as pd


# alphabet, global variable
alphabet = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']


def create_msa_array(msa):
    '''
        remove columns where the subject sequence has gaps, the top sequence is assumed
        to be the subject sequence and gaps are assumed to be represented as '-'s
    '''
    print('Computing a NumPy ndarray representation of the alignment')
    msa_arr = np.asarray([list(rec.upper()) for rec in msa])
    msa_arr[msa_arr == 'B'] = 'D'
    msa_arr[msa_arr == 'Z'] = 'Q'
    return msa_arr


def compute_aa_probs(column_i):
    '''
    '''
    # amino acid background relative frequencies 
    # bg_rf = np.array([0.085786, 0.045676, 0.047306, 0.058022, 0.018036, 0.037722, 0.059724, 0.081155, 0.021639, 
    #        0.052944, 0.081156, 0.058717, 0.021109, 0.039946, 0.048178,0.063047, 0.060835, 0.014256, 0.036310, 0.068436])
    f = np.asarray([float(sum(column_i == aa) + 1) for aa in alphabet])
    rf = np.true_divide(f, sum(f))
    return rf


def compute_joint_probs(column_i, column_j):
    '''
    Compute the joint distribution of amino acid pairs for position pair i, j.
    See Weigt et al, PNAS, 2009 for technical details.
    '''
    if len(column_i) != len(column_j):
        print('Unequal depth of the two given alignment column pair.')
        sys.exit(0)
    
    # compute joint probabilities
    joint_probs = []
    for a in alphabet:
        # count frequency of amino acid pair (a, b) in columns (i, j)
        joint_freq = np.asarray([sum((column_i == a) * (column_j == b)) + (1. / 20) for b in alphabet])

        # compute joint probability
        joint_prob = joint_freq / (20 + len(column_i))
        joint_probs.append(joint_prob)

    return np.asarray(joint_probs)


def safe_log(x):
    '''
    '''
    return 0 if x <= 0 else math.log(x)


def compute_mi(probs_i, probs_j, joint_probs_ij):
    '''
    Compute the mutual information between position pair(i, j) from the given multiple sequence alignment.
    See Weigt et al, PNAS, 2009 for technical details.
    '''
    mi = 0
    for a in range(20):
        for b in range(20):
            mi += joint_probs_ij[a, b] * safe_log(joint_probs_ij[a, b] / (probs_i[a] * probs_j[b]))
    return mi


def compute_excess_mi(column_i, column_j):
    '''
    Correction of mutual information for sampling bias.
    See Cline et al, Proteins, 2002.
    '''
    # discard pairs within a given sequence where either i or j is a gap
    valid_pairs = (column_i != '-') * (column_j != '-')
    column_i = column_i[valid_pairs]
    column_j = column_j[valid_pairs]

    # compute amino acid probabilities at each column
    probs_i = compute_aa_probs(column_i)
    probs_j = compute_aa_probs(column_j)

    joint_probs = compute_joint_probs(column_i, column_j)
    apparent_mi = compute_mi(probs_i, probs_j, joint_probs)
    print('The apparent mutual information is:', apparent_mi)

    # permutation test for significance
    n_permutations = 1000
    empirical_mis = []
    for i in range(n_permutations):
        column_j_permuted = np.random.permutation(column_j)
        joint_probs = compute_joint_probs(column_i, column_j_permuted)
        mi_permutation = compute_mi(probs_i, probs_j, joint_probs)
        print('Mutual information from permutation round', i, 'is:', mi_permutation)
        empirical_mis.append(mi_permutation)
    n = sum(np.asarray(empirical_mis) > apparent_mi)
    if n >= 10:
        print('Got', n, 'permutations where the mutual information is greater than the apparent mutual information:', apparent_mi, 'failed to reject the null hypothesis.')
        return 0.0
    else:
        # correct for bias
        return apparent_mi - np.mean(empirical_mis)


def compute_shortest_distance(res_a, res_b):
    '''
    Returns the smallest of the distances between the heavy atoms of residue a and those of residue b.
    '''
    distances = [np.linalg.norm(atom_a - atom_b) for atom_a in res_a.get_list() for atom_b in res_b.get_list()]
    return min(distances)


def main():
    '''
    '''
    # parser for parsing command line arguments
    parser = ArgumentParser(description='Compute mutual information from a multiple sequence alignment for a given list of pairs of residues.')
    parser.add_argument('-a', '--msa', dest='msa', help='multiple sequence alignment in fasta format')
    parser.add_argument('-p', '--pdb', dest='pdb', help='input file in PDB format')
    parser.add_argument('-r', '--res', dest='res', nargs=2, help='two files, each containing a list of residue numbers')
    parser.add_argument('-d', '--dm', dest='dm', help='a matrix of pairwise distances')
    parser.add_argument('-m', '--mi', dest='mi', help='name to the file to which to write the mutual information')

    # parse command line arguments
    args = parser.parse_args()

    # read the multiple sequence alignment
    alignment = AlignIO.read(args.msa, format='fasta')

    # create a 2D numpy array from the alignment
    msa_array = create_msa_array(alignment)

    # parse the given PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='pdb', file=args.pdb)
    model = structure[0]

    # residues given on the command line
    with open(args.res[0], 'rt') as f:
        # residue IDs contain both chain ID and sequence position information
        res_ids_a = [tuple(re.split(r'[,;\s]\s*', line.strip())) for line in f.readlines()]
    res_pos_a = [int(i[1]) for i in res_ids_a]
    chain_a = res_ids_a[0][0]
    residues_a = [model[chain_a][int(i[1])] for i in res_ids_a]

    with open(args.res[1], 'rt') as f:
        # residue IDs contain both chain ID and sequence position information
        res_ids_b = [tuple(re.split(r'[,;\s]\s*', line.strip())) for line in f.readlines()]
    res_pos_b = [int(i[1]) for i in res_ids_b]
    chain_b = res_ids_b[0][0]
    residues_b = [model[chain_b][int(i[1])] for i in res_ids_b]

    # compute residue pair distances
    pair_distances = np.asarray([[compute_shortest_distance(a, b) for b in residues_b] for a in residues_a])
    
    # compute the mutual information between the distributions of amino acids occurring in each pair of positions
    computed_pairs = {}
    mi = []
    for i, a in enumerate(res_pos_a):
        mi_a = []
        column_a = msa_array[:, a-1].copy()
        gap_fraction_a = sum(column_a == '-') / len(column_a)
        # ignore columns where there are more than 25% gaps
        if gap_fraction_a > 0.25:
            mi_a = [np.NAN for i in range(len(residues_b))]
            mi.append(mi_a)
            continue
        for j, b in enumerate(res_pos_b):
            column_b = msa_array[:, b-1].copy()
            gap_fraction_b = sum(column_b == '-') / len(column_b)
            # ignore columns where there are more than 25% gaps
            if gap_fraction_b > 0.25 or abs(a - b) <= 4: # or pair_distances[i][j] > 4.5:
                mi_a.append(np.NAN)
                continue
            # mutual information is symmetric: I(X, Y) = I(Y, X)
            if (b, a) in computed_pairs.keys():
                print('The mutual information between', (a, b), 'was already computed.')
                x, y = computed_pairs[(b, a)]
                mi_a.append(mi[x][y])
                continue
            print('Now computing the mutual information between column pair: ' + '(' + str(a) + ', ' + str(b) + ')')
            mi_a.append(compute_excess_mi(column_a, column_b))
            computed_pairs[(a, b)] = (i, j)
        mi.append(mi_a)

    # make a DataFrame for formatting output
    subject_seq = msa_array[0]
    indices = [''.join([i[0], i[1], subject_seq[int(i[1])-1]]) for i in res_ids_a]
    headers = [''.join([i[0], i[1], subject_seq[int(i[1])-1]]) for i in res_ids_b]
    mi_df = pd.DataFrame(np.array(mi), index=indices, columns=headers)
    pair_distances_df = pd.DataFrame(pair_distances, index=indices, columns=headers)

    # write to file
    pair_distances_df.to_csv(args.dm, index=True, header=True, sep=',', float_format='%.2f')
    mi_df.to_csv(args.mi, index=True, header=True, sep=',', float_format='%.4f')


if __name__ == '__main__':
    main()
