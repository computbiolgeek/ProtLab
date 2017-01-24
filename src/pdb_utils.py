"""
    @summary: 
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu or comput.biol.geek@gmail.com
    @change: 
    @copyright: 
"""

import numpy as np
from collections import defaultdict

def ReadPDB(pdb_file, atom="CA"):
    """
        @summary: Takes as input a PDB file, extracts residue names and the
            corresponding alpha carbon coordinates.
        @param pdb_file: the string filename of a PDB file
        @return: 
            residue names: a length-N list containing the three-letter code for
                the type of each of the corresponding amino acid residues
            coordinates: a (N, 3) dimensional NumPy array containing the position
                of the alpha carbon in each of the N amino acid residues         
    """
    res_names = []
    coordinates = []
    with open(pdb_file, "rt") as f:
        for line in f:
            if line.startswith("ATOM"):
                if atom in line[12:16]:
                    res_names.append(line[17:20])
                    coordinates.append([float(line[30:38]),
                                              float(line[38:46]),
                                              float(line[46:54])])
            elif line.startswith("TER"):
                break
    return res_names, np.array(coordinates)

def ResHydrophobic(residues):
    """
        @summary: Takes a length-N list of three-letter residue codes, returns
            a length-N array with Boolean values indicating whether the corresponding
            residue is hydrophobic or not.
        @param residues: a length-N list of three-letter residue codes
        @return:
            is_hydrophobic, a length-N NumPy array with True or False boolean
            values indicating whether the residue codes in "residues" correspond
            hydrophobic or hydrophilic residues
    """
    hydrophobic = {"ALA", 
                   "CYS", 
                   "PHE", 
                   "ILE", 
                   "LEU", 
                   "MET", 
                   "VAL"}
    is_hydrophobic = [residue in hydrophobic for residue in residues]
    return np.array(is_hydrophobic)
                    
def RadiusOfGyration(coordinates):
    """
        @summary: Computes the radius of gyration of the passed set of atoms.
        @param 
            coordinates: a (N, 3) dimensional NumPy array.
        @return: 
            Rg: a float with the corresponding radius of gyration 
    """ 
    r_mean = coordinates.mean(axis=0)
    rg = np.mean(np.sum((coordinates - r_mean)**2))
    return np.sqrt(rg)  

def GetContacts(coordinates, residues, cutoff):
    """
        @summary: 
        @param 
            coordinates:
            residues:
            cutoff:
        @return:
            res_freqs:
            contacts:
    """
    amino_acids = ["ALA", 
                   "CYS",
                   "ASP",
                   "GLU", 
                   "PHE",
                   "GLY",
                   "HIS", 
                   "ILE",
                   "LYS",
                   "LEU",
                   "MET",
                   "ASN",
                   "PRO",
                   "GLN",
                   "ARG",
                   "SER",
                   "THR",
                   "VAL",
                   "TRP",
                   "TYR"]
    res_freqs = {res: residues.count(res) for res in amino_acids}
    contacts = defaultdict(list)
    # number of residues
    n = coordinates.shape[0]
    # compute the upper triangular matrix
    for i in range(n):
        for j in range(i + 1, n):
            dist_ij = np.sqrt((sum(coordinates[i] - coordinates[j])**2))
            if dist_ij < cutoff:
                contacts[(residues[i], residues[j])].append(dist_ij)
    return res_freqs, contacts