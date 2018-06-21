#!/home/bian/anaconda3/bin/python3

import os, re
from collections import defaultdict, OrderedDict
import numpy as np
import csv
from Bio.PDB import PDBParser, NeighborSearch
from argparse import ArgumentParser


RESIDUE_POLARITY_TYPE = {
    'PHE': 0,
    'MET': 0,
    'ILE': 0,
    'LEU': 0,
    'VAL': 0,
    'CYS': 0,
    'PRO': 0,
    'ALA': 0,
    'GLY': 0,
    'THR': 0.5,
    'TRP': 0.5,
    'SER': 0.5,
    'TYR': 0.5,
    'GLN': 0.5,
    'ASN': 0.5,
    'HIS': 1,
    'GLU': -1,
    'LYS': 1,
    'ASP': -1,
    'ARG': 1 
}

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


class Contact:
    '''
    '''
    def __init__(self, res_a, res_b):
        '''
        Parameters
        ----------
        res_a : Residue
            The first residue making this contact. A Biopython Residue object.
        res_b : Residue
            The second residue making this contact. A Biopython Residue object.
        '''
        self._res_a = res_a
        self._res_b = res_b

    def __str__(self):
        '''
        Returns
        -------
        str
            A user friendly representation of this contact.
        '''
        return self._res_a.get_resname() + str(self._res_a.get_id()[1]) + '--' \
            + self._res_b.get_resname() + str(self._res_b.get_id()[1])

    def get_res_a(self):
        '''
        Returns
        -------
        Residue
            The first residue of this contact.
        '''
        return self._res_a

    def get_res_b(self):
        '''
        Returns
        -------
        Residue
            The second residue of this contact.
        '''
        return self._res_b

    def get_res_a_index(self):
        '''
        Returns
        -------
        int
            Sequence index of the first residue.
        -------

        '''
        return self._res_a.get_id()[1]

    def get_res_b_index(self):
        '''
        Returns
        -------
        int
            Sequence index of the second residue.
        '''
        return self._res_b.get_id()[1]

    def get_shortest_distance(self):
        '''
        Returns
        -------
        float
            The smallest of the distances between the heavy atoms
            of residue a and those of residue b.
        '''
        distances = [atom_a - atom_b for atom_a in self._res_a.get_list() 
                for atom_b in self._res_b.get_list()]
        return min(distances)

    def get_cb_distance(self):
        '''
        Returns
        -------
        float
            The distance between the Cb atoms of the contacting residues.
        '''
        try:
            res_a_cb = self._res_a['CB']
        except KeyError:
            res_a_cb = self._res_a['CA']
        try:
            res_b_cb = self._res_b['CB']
        except KeyError:
            res_b_cb = self._res_b['CA']
        return res_a_cb - res_b_cb

    @staticmethod
    def get_sidechain_centroid(res):
        '''
            C = (v1 + ... + vn) / n
        '''
        # return CA coordinates if residue is GLY
        if res.get_resname() == 'GLY':
            return res['CA'].get_coord()
        # for other residues, return sidechain centroid
        centroid = np.array( [0.0, 0.0, 0.0] )
        number_sidechain_atoms = 0
        for atom in res.get_list():
            if atom.get_name() not in BACKBONE_ATOMS:
                centroid += atom.get_coord()
                number_sidechain_atoms += 1
        return centroid / number_sidechain_atoms

    def get_centroid_distance(self):
        '''
        Returns
        -------
        float
            The distance between the side-chain centroids of the contacting residues.
        '''
        centroid_a = self.get_sidechain_centroid(self._res_a)
        centroid_b = self.get_sidechain_centroid(self._res_b)
        return np.linal.norm(centroid_a - centroid_b)

    def is_electrostatic(self):
        '''
        Determines whether the contact is of is_electrostatic nature.

        Returns
        -------
        bool
            True if the contacting residues are oppositely charged.
        '''
        res_a_type = RESIDUE_POLARITY_TYPE[self._res_a.get_resname()]
        res_b_type = RESIDUE_POLARITY_TYPE[self._res_b.get_resname()]
        return res_a_type * res_b_type < 0

    def is_nonpolar(self):
        '''
        Determines whether the contact is of nonpolar nature.

        Returns
        -------
        bool
            True if none of the contact residues is polar.
        '''
        res_a_type = RESIDUE_POLARITY_TYPE[self._res_a.get_resname()]
        res_b_type = RESIDUE_POLARITY_TYPE[self._res_b.get_resname()]
        return res_a_type == 0 and res_b_type == 0

    def is_polar(self):
        '''
        Determines whether the contact is of polar nature.
        
        Returns
        -------
        bool
            True if the either of the contact residues is polar.
        '''
        return not self.is_nonpolar() 

    def is_local(self):
        '''
        @TODO
        '''
        pass

    def is_global(self):
        '''
        Determines whether a contact is global.

        Returns
        -------
        bool
            True if the contact is global.
        '''
        return not self.is_local()


def get_residues(chain, start, end):
    '''
    '''
    residues = []
    for i in range(start, end+1):
        try:
            residues.append(chain[i])
        except KeyError:
            continue
    return residues


def search_for_all_contacts(residues, radius=4.5):
    '''
    Search for all contacts in the given set of residues based on 
    distances between heavy atoms (atoms other than hydrogens).
    '''
    atom_list = []
    for r in residues:
        if r.get_resname() == 'GLY':
            atom_list.append(r['CA'])
        else:
            # atom_list += [a for a in r.get_atoms() if a.get_name() not in BACKBONE_ATOMS]
            atom_list += [a for a in r.get_atoms() if a.get_name() == 'CB']
    ns = NeighborSearch(atom_list)
    all_contacts = [Contact(res_a = c[0], res_b = c[1]) for c in ns.search_all(radius, level='R')]
    return all_contacts


def count_local_global_contacts(contacts, idx_to_seg_id_mapping):
    '''
    Computes the number of local contacts and the number of global contacts.
    A contact is local if it is made by residues from sequentially neighboring
    segments, and global if it is made by residues from sequentially separated
    segments. Ignore contacts made by residues in the same segment.

    Parameters
    ----------
    contacts : list 
        A list of contacts identified through calling the search_all() method
        on a NeighborSearch object.
    idx_to_seg_id_mapping : dict
        A dictionary containing the mapping from residue indices from segment IDs.

    Returns
    -------
    tuple
        Number of local contacts and number of global contacts
    '''
    num_local = 0
    num_global = 0
    for c in contacts:
        res_id_a = c.get_res_a_index()
        seg_id_a = idx_to_seg_id_mapping[res_id_a]
        res_id_b = c.get_res_b_index()
        seg_id_b = idx_to_seg_id_mapping[res_id_b]
        # local contacts if seg_id_a and seg_id_b differs by 1
        # ignore contacts made by residues in the same segment
        if abs(seg_id_a - seg_id_b) == 1:
            num_local += 1
        elif abs(seg_id_a - seg_id_b) >= 2:
            num_global += 1
    return num_local, num_global


def map_contacts_to_segments(contacts, idx_to_seg_id_mapping):
    '''
    Map contacts to segments based on the indices of residues making the contact.

    Parameters
    ----------
    contacts : list
        A list of contacts identified by the Biopython NeighborSearch method. Each contact
        is a tuple of a pair of residues making the contact.
    idx_to_seg_id_mapping: dict
        A mapping of residue indices to segment IDs.

    Returns
    -------
    dict
        A mapping of contacts to segment IDs.
    '''
    contacts_to_segments_mapping = {}
    for c in contacts:
        res_id_a = c.get_res_a_index()
        seg_id_a = idx_to_seg_id_mapping[res_id_a]
        res_id_b = c.get_res_b_index()
        seg_id_b = idx_to_seg_id_mapping[res_id_b]
        # exclude contacts made by residues from within the same segment
        if seg_id_a != seg_id_b:
             contacts_to_segments_mapping[c] = (seg_id_a, seg_id_b)
    return contacts_to_segments_mapping


def summarize_segment_contacts(contacts_to_segments_mapping):
    '''
    Count the total number of contacts, the number of local, global, nonpolar,
    polar, electrostatic contacts, respectively, of each segment.

    Parameters
    ----------
    contacts_to_segments_mapping : dict
        A dictionary that maps contacts to segment pairs. The key of the dictionary is
        residue pair (contact) and the value is the corresponding segment pair.

    Returns
    -------
    OrderedDict
        An odered dictionary that maps segment IDs to summary of contacts. The key of
        the dictionary is segment ID and the value is the corresponding summary of 
        contacts stored in an ordered dictionary.
    '''
    # get all contacts for each segment
    segment_contacts = defaultdict(list)
    for k, v in contacts_to_segments_mapping.items():
        segment_contacts[v[0]].append(k)
        segment_contacts[v[1]].append(k)

    # summarize the contacts for each segment
    segment_contacts_summary = OrderedDict() 
    for seg_id, contacts in sorted(segment_contacts.items()):
        total_contacts = len(contacts)
        global_contacts = 0
        local_contacts = 0
        polar_contacts = 0
        nonpolar_contacts = 0
        global_polar_contacts = 0
        local_polar_contacts = 0
        electrostatic_contacts = 0
        for c in contacts:
            seg_pair = contacts_to_segments_mapping[c]
            if abs(seg_pair[0] - seg_pair[1]) >= 2:
                global_contacts += 1
            else:
                local_contacts += 1
            if c.is_nonpolar():
                nonpolar_contacts += 1
            else:
                polar_contacts += 1
                if c.is_electrostatic():
                    electrostatic_contacts += 1
            if abs(seg_pair[0] - seg_pair[1]) >= 2 and not c.is_nonpolar():
                global_polar_contacts += 1
            if abs(seg_pair[0] - seg_pair[1]) == 1 and not c.is_nonpolar():
                local_polar_contacts += 1
        contacts_summary = [
                total_contacts,
                global_contacts,
                local_contacts,
                polar_contacts,
                nonpolar_contacts,
                global_polar_contacts,
                global_contacts - global_polar_contacts,
                local_polar_contacts,
                local_contacts - local_polar_contacts,
                electrostatic_contacts]

        segment_contacts_summary[seg_id] = contacts_summary
    return segment_contacts_summary


def parse_segments(segment_file):
    '''
    Parse the segment file and return a list of tuples whose three elements 
    are segment ID and residue index ranges for each segment.
    '''
    id_and_idx_ranges = []
    with open(segment_file, 'rt') as f:
        for l in f:
            segment = map(int, re.split(r'[:-]', l.strip()))
            id_and_idx_ranges.append(tuple(segment))
    return id_and_idx_ranges


def get_model_from_pdb(pdb_file, model_number=0):
    '''
    Parse the PDB file and return a given structure model.
    '''
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(id=os.path.basename(pdb_file).split('.')[0],
            file=pdb_file)
    model = structure[model_number]
    return model

def map_index_to_segment_id(segment_tuples):
    '''
    Return a dict containing the mapping from residue indices to segment IDs.
    '''
    idx_to_seg_id = {}
    for s in segment_tuples:
        for i in range(s[1], s[2]+1):
            idx_to_seg_id[i] = s[0]
    return idx_to_seg_id


def parse_commandline_arguments():
    '''
    Parse command line arguments and return them.
    '''
    arg_parser = ArgumentParser()
    arg_parser.add_argument('-p', '--pdb', dest='pdb', required=True, help='input PDB file')
    arg_parser.add_argument('-s', '--segments', dest='segments', required=True, help='file containing segment IDs and residue index ranges')
    arg_parser.add_argument('-r', '--radius', dest='radius', type=float, default=4.5, help='radius to use when searching for contacts')
    arg_parser.add_argument('-c', '--contacts', dest='contacts', required=True, help='file to which to write residue contacts')
    arg_parser.add_argument('-o', '--output', dest='output', help='file to which to write a summary of contacts')
    return arg_parser.parse_args()


def main():
    '''
    '''
    # parse command line arguments
    args = parse_commandline_arguments()    

    # parse the PDB file to get a structure model
    model = get_model_from_pdb(args.pdb)

    # consider only the fist chain
    chain = next(model.get_chains())

    # parse the segment file
    segment_tuples = parse_segments(args.segments)

    # create a dict that maps residue indices to helix IDs
    idx_to_helix_id = map_index_to_segment_id(segment_tuples)

    # get all interested residues
    residues = []
    for s in segment_tuples:
        residues += get_residues(chain, s[1], s[2])

    # search for all contacts
    all_contacts = search_for_all_contacts(residues, radius=args.radius)
    
    # map contacts to segments
    contacts_to_segments_mapping = map_contacts_to_segments(all_contacts, idx_to_helix_id)
    
    # get all contacts for each segment
    segment_contacts = defaultdict(list)
    for k, v in contacts_to_segments_mapping.items():
        segment_contacts[v[0]].append(k)
        segment_contacts[v[1]].append(k)

    # get the PDB id
    pdb_id = os.path.basename(args.pdb).split('.')[0]
    
    # write segment contacts to file
    with open(args.contacts, 'wt') as f:
        for i, contacts in sorted(segment_contacts.items()):
            f.write(pdb_id + ',' + str(i) + ',' + ','.join(str(c) for c in contacts) + '\n')

    # summarize contacts for each segment
    segment_contacts_summary = summarize_segment_contacts(contacts_to_segments_mapping)

    # write summary to file
    with open(args.output, 'wt') as f:
        csv_writer = csv.writer(f)
        headers = ['pdb id', 'segment id', 'segment_range', 'total contacts', 
                'global contacts', 'local contacts','polar contacts', 
                'nonpolar contacts', 'global polar contacts', 'global nonpolar contacts',
                'local polar contacts', 'local nonpolar contacts',
                'electrostatic contacts']
        csv_writer.writerow(headers)
        seg_ids_contacts = list(segment_contacts_summary.keys())
        for seg_id in range(1, segment_tuples[-1][0] + 1):
            if seg_id not in seg_ids_contacts:
                csv_writer.writerow([pdb_id, seg_id] + ['%s-%s' % segment_tuples[seg_id - 1][1:]] + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            else:
                csv_writer.writerow([pdb_id, seg_id] + ['%s-%s' % segment_tuples[seg_id - 1][1:]] + segment_contacts_summary[seg_id])



if __name__ == '__main__':
    main()
