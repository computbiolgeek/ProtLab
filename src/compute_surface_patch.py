#!/usr/bin/env python3

import numpy as np
import sys
import re
import math
import heapq
from argparse import ArgumentParser
from Bio.PDB import PDBParser, NeighborSearch


class SurfacePatch:
    '''
    '''
    # kyte-doolittle hydrophobicity scale
    kyte_hp_scale = {
        'PHE': 2.8,
        'MET': 1.9,
        'ILE': 4.5,
        'LEU': 3.8,
        'VAL': 4.2,
        'CYS': 2.5,
        'ALA': 1.8,
        'THR': -0.7,
        'TRP': -0.9,
        'GLY': -0.4,
        'SER': -0.8,
        'PRO': -1.6,
        'TYR': -1.3,
        'HIS': -3.2,
        'GLN': -3.5,
        'ASN': -3.5,
        'GLU': -3.5,
        'LYS': -3.9,
        'ASP': -3.5,
        'ARG': -4.5 
    }

    def __init__(self, residues):
        self.residues = residues

    def get_mean_hydrophobicity(self):
        '''
        Calculates the mean hydrophobicity of the patch.
        '''
        hp = 0
        for res in self.residues:
            hp += self.kyte_hp_scale[res.get_resname()]
        return hp / len(self.residues)
    
    def get_solvation_potential(self):
        '''
        '''
        pass

    def get_mean_rsa(self):
        '''
        Calculates the mean relatie accessible surface area of the patch.
        '''
        rsa = 0
        for res in self.residues:
            rsa += res.rsa
        return rsa / len(self.residues)
    
    def get_planarity(self):
        '''
        '''
        pass



class SurfacePatchScanner:
    '''
    '''

    def __init__(self, model, rsas, size):
        '''
        model : an object of the Biopython Model class
        size : int
        rsas : dict
        '''
        self.model = model
        self.size = size
        self.rsas = rsas
        self.patches = []

    def scan_for_patches(self):
        '''
        '''
        ca_list = [res['CA'] for res in self.get_surface_residues() if 'CA' in res]
        neighbors_distances = dict()
        sol_vecs = dict()
        neighbor_search = NeighborSearch(ca_list)
        for ca in ca_list:
            # search for neighbors of each residue and store them 
            neighbors = neighbor_search.search(ca.coord, math.inf, level='A')
            distances = neighbor_search.kdt.get_radii()
            current_neighbors_distances = dict(zip(neighbors, distances))
            neighbors_distances[ca] = current_neighbors_distances
            
            # compute the solvent vector for each surface residue
            ten_nn = heapq.nsmallest(10, current_neighbors_distances, 
                    key=lambda k: current_neighbors_distances[k])
            sol_vecs[ca] = ca.coord - np.mean([a.coord for a in ten_nn], axis=0)
        
        already_scanned = set()
        for ca in ca_list:
            if ca in already_scanned:
                continue
            
            # get k nearest neighbors where k = size
            knn = heapq.nsmallest(self.size, neighbors_distances[ca], key=lambda k: neighbors_distances[ca][k])
            
            # remove neighbors that do not satisfy requirements
            knn = [n for n in knn if self._satisfy_requirement(sol_vecs[ca], sol_vecs[n])]
            
            # add residues into the list of patches
            self.patches.append(SurfacePatch([a.get_parent() for a in knn]))
            already_scanned.update(knn)

    def get_num_patches(self):
        '''
        '''
        return len(self.patches)

    def get_surface_residues(self):
        '''
        '''
        for res in self.model.get_residues():
            res_pos = res.get_id()[1]
            if self.rsas[res_pos] > 0.05:
                yield res
    
    def _satisfy_requirement(self, vec_a, vec_b):
        '''
        '''
        return np.dot(vec_a, vec_b) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_b)) > -0.3420

    def __iter__(self):
        '''
        '''
        return iter(self.patches)

    def __getitem__(self, key):
        '''
        '''
        return self.patches[key]
    

def main():
    '''
    '''
    # create an object for parsing commandline arguments
    arg_parser = ArgumentParser(description='This script computes protein surface patches given a PDB file and a file of residue accessibility.'
        'See Jones and Thornton, JMB, 1997 for technical details.'
    )
    arg_parser.add_argument('-p', '--pdb', dest='pdb', required=True, help='PDB file containing coordinates of protein atoms')
    arg_parser.add_argument('-a', '--rsa', dest='rsa', required=True, help='file containing residue relative solvent accessibilities')
    arg_parser.add_argument('-n', '--size', dest='size', type=int, required=True, help='the size (number of residues) of the surface patches')

    # parse commandline arguments
    args = arg_parser.parse_args()

    # get the structure/model
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(id='prot', file=args.pdb)
    model = structure[0] # use the first model in the PDB file

    # get residue accessible surface areas
    rsas = {}
    with open(args.rsa, 'rt') as f:
        for line in f:
            if not line.strip().startswith('#'):
                seq_pos, rsa = re.split(r'[,:\s]\s*', line.strip())
                rsas[int(seq_pos)] = float(rsa) / 100

    # add a rsa attribute to each residue
    for res in model.get_residues():
        res.rsa = rsas[res.get_id()[1]]
    
    patch_scanner = SurfacePatchScanner(model, rsas, args.size)
    patch_scanner.scan_for_patches()
    print('Number of patches:', patch_scanner.get_num_patches())

    for patch in patch_scanner:
        print(patch.get_mean_rsa())
    
    
if __name__ == '__main__':
    main()
