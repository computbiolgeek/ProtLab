import numpy as np

# read bit-valued substitution matrices into dictionaries
filename = '/hd0/lib14/workspace/seq_utils/data/PAM250.txt'
with open(filename, 'rt') as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines if not line.startswith('#')]

aas = lines.pop(0)
aas = aas.split()
ll = [[int(i) for i in line.split()[1:]] for line in lines]
# np.array(ll)
dd = {line.split()[0]: dict(zip(aas, [int(i) for i in line.split()[1:]])) for line in lines}
# write dd into file
with open('pam250.txt', 'wt') as f:
    f.write('{')
    items = list(dd.items())
    [f.write("'" + str(key) + "': " + str(value) + ",\n") for key,value in items[:-1]]
    f.write("'" + str(key) + "': " + str(value) + "}")
# od = OrderedDict(sorted(dd.items()))
# arr = np.array([OrderedDict(sorted(d.items())).values() for d in od.values()])


# read BLOSUM clustered target frequencies qij in to a symmetric matrix
filename = '/hd0/lib14/workspace/seq_utils/data/blosum50.qij'
with open(filename, 'rt') as f:
    lines = f.readlines()
    line = ''
    for l in lines:
        if not l.startswith('#') and l.find('A') == -1:
            line += l

freqs = [float(i) for i in line.split()]
m = np.zeros((20,20))
m[np.tril_indices(20,0)] = freqs
m += m.T
m[np.diag_indices(20)] = np.diag(m)/2

# def read_fasta_alignment(filename):
#     """ Read in the alignment stored in the FASTA file, filename. Return two
#     lists: the identifiers and sequences. """
# 
#     f = open(filename)
# 
#     names = []
#     alignment = []
#     cur_seq = ''
# 
#     for line in f:
#         line = line[:-1]
#         if len(line) == 0:
#             continue
# 
#         if line[0] == ';': 
#             continue
#         if line[0] == '>':
#             names.append(line[1:].replace('\r', ''))
#             if cur_seq != '':
#                 cur_seq = cur_seq.upper()
#                 for i, aa in enumerate(cur_seq):
#                     if aa not in iupac_alphabet:
#                         cur_seq = cur_seq.replace(aa, '-')
#                 alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
#                 cur_seq = ''
#         elif line[0] in iupac_alphabet:
#             cur_seq += line.replace('\r', '')
# 
#     # add the last sequence
#     cur_seq = cur_seq.upper()
#     for i, aa in enumerate(cur_seq):
#         if aa not in iupac_alphabet:
#             cur_seq = cur_seq.replace(aa, '-')
#     alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
# 
#     return names, alignment


from substitution_matrices import BLOSUM50_QIJ
BLOSUM50_QIJ