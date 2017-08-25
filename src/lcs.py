"""
    @summary: Longest common subsequence problem. Find the longest subsequence common
        to two strings. A subsequence of a string is simply an (ordered) sequence of
        characters (not necessarily consecutive).
    @param param: Two strings, v and w.
    @return: The longest common subsequence of v and w.
"""

def lcs(v, w):
    """
    """
    len_v = len(v)
    len_w = len(w)
    # s is the dynamic programming table
    s = [[None]*(len_w + 1)]*(len_v + 1)
    for i in range(len_v + 1):
        for j in range(len_w + 1):
            if i == 0 or j == 0:
                s[i][j] = 0
            elif v[i - 1] ==  w[j - 1]:
                s[i][j] = s[i - 1][j - 1] + 1
            else:
                s[i][j] = max(s[i][j - 1], s[i - 1][j])
    return s[len_v][len_w]

def needleman_wunsch(v, w, m):
    """
        Summary:
            Find the best alignment between two sequences under a given scoring matrix.
        
        Parameter(s):
            v(str): The first sequence.
            w(str): The second sequence.
            m(dict): The scoring matrix. 
        
        Return:
            An alignment of v and w whose score is maximal among all possible alignments
            of v and w under the given matrix m.
        
        See also:
            Needleman and Wunsch, JMB, 1970
    """
    len_v = len(v)
    len_w = len(w)
    # s is the dynamic programming table
    s = [[None]*(len_w + 1)]*(len_v + 1)
    # gap penalty, typically a negative number in a substitution matrix
    gp = m['A']['*']
    # fill the first row
    for i in range(len_w + 1):
        s[0][i] = gp * i
    # fill the first column
    for j in range(1, len_v + 1):
        s[j][0] = gp * j
    # fill the rest of the dynamic programming table
    for i in range(1, len_v + 1):
        for j in range(1, len_w + 1):
            s[i][j] = max(s[i - 1][j] + gp, s[i][j - 1] + gp, s[i - 1][j - 1] + m[v[i - 1]][w[j - 1]])
    # return the maximum match
    return s[len_v][len_w]


def read_fasta(filename):
    """
    """
    ids = []
    seqs = []
    with open(filename, 'rt') as f:
        cur_seq = ''
        for line in f.readlines():
            if line.startswith('>'):
                ids.append(line.strip())
                if cur_seq != '':
                    cur_seq = cur_seq.upper()
                    seqs.append(cur_seq)
                    cur_seq = ''
            else:
                cur_seq += line.strip()
        # add the last sequence
        cur_seq = cur_seq.upper()
        seqs.append(cur_seq)
    return ids, seqs