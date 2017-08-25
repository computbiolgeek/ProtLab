"""
    @summary: 
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: 
    @copyright: 
"""

import numpy as np
        
class MultipleAlignmentArray( object ):
    """
    """
    
    alphabet = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    # amino acid background relative frequencies 
    bg_rf = np.array([
                        0.085786, 
                        0.045676, 
                        0.047306, 
                        0.058022, 
                        0.018036,
                        0.037722, 
                        0.059724, 
                        0.081155, 
                        0.021639, 
                        0.052944,
                        0.081156, 
                        0.058717, 
                        0.021109, 
                        0.039946, 
                        0.048178,
                        0.063047, 
                        0.060835, 
                        0.014256, 
                        0.036310, 
                        0.068436
                    ])
    
    groups = [
        ['A','V','L','I','M','C'], # aliphatic 
        ['F','W','Y','H'],  # aromatic 
        ['S','T','N','Q'],  # polar 
        ['K','R'],  # positive 
        ['D','E'],  # negative 
        ['G','P'] # glycine and proline
        ]    
        
    def __init__( self, msa = None, subject_index = 0, gap = '-', pseudocount = 0 ):
        self.__msa = msa
        self.__subject_index = subject_index
        self.__gap = gap
        self.__msa_arr = None
        self.__ncols = self.msa_arr.shape[1]
        self.__nrows = self.msa_arr.shape[0]
        self.__pc = pseudocount
        self.__f = None # amino acid frequencies 
        self.__rf = None 
        self.__wf = None # weighted relative frequencies
        self.__icf = None
        self.__gap_f = None
        self.__gap_rf = None
        self.__group_f = None
        self.__group_rf = None
        self.__aa_in_col = None 
        self.__pssm = None
        
    @property
    def subject_index( self ):
        """
        """
        return self.__subject_index
    
    @subject_index.setter
    def subject_index( self, index ):
        """
        """
        if index < 0 or index > len( self.__msa ):
            raise ValueError( "Invalid index of subject sequence" )
        self.__subject_index = index
        
    @property
    def gap( self ):
        """
        """
        return self.__gap
    
    @gap.setter
    def gap( self, gap ):
        """
        """
        if gap not in ['-', '.', '*']:
            raise ValueError( "Invalid symbol for gap" )
        self.__gap = gap
           
    @property
    def msa_arr( self ):
        """
            remove columns where the subject sequence has gaps, the top sequence is assumed
            to be the subject sequence and gaps are assumed to be represented as '-'s
        """
        if self.__msa_arr is None:
            print( "Computing a NumPy ndarray representation of the alignment" )
            self.__msa_arr = np.array( [list( rec.upper() ) for rec in self.__msa] )
            print( "Removing gaps in the subject sequence" )
            # set X to gap
            self.__msa_arr[self.__msa_arr == 'X'] = self.gap
            self.__msa_arr[self.__msa_arr == 'B'] = 'D'
            self.__msa_arr[self.__msa_arr == 'Z'] = 'Q'
            self.__msa_arr = self.__msa_arr[:, self.__msa_arr[self.__subject_index] != self.__gap ]
        return self.__msa_arr
    
    @property
    def f( self ):
        """

        """
        if self.__f is None:
            print( "Counting occurrence of each amino acid type in each aligned column." )
            self.__f = np.array( [[float(sum( self.msa_arr[:, i] == aa )) for aa in MultipleAlignmentArray.alphabet] 
                         for i in range( self.__ncols )] )
            # add pseudocount to the observed frequency of each amino acid in each aligned column
            self.__f += (self.__pc * MultipleAlignmentArray.bg_rf) 
        return self.__f
    
    @property
    def rf( self ):
        """
            @note: amino acid frequencies relative to the total number of 
            non-gapped sequences at each position
        """
        if self.__rf is None:
            print( "Calculating the probability that each amino acid type in each aligned column." ) 
            self.__rf = np.true_divide( self.f, self.f.sum( axis = 1 )[:, np.newaxis] )          
        return self.__rf
    
    @property
    def wf(self):
        """
            @summary: weighted frequencies
            @see: Pei, 2001, Bioinformatics
        """
        if self.__wf is None:
            all_wf = []
            sw = self.sequence_weights_henikoff()
            # iterate over aligned columns
            for i in range(self.__ncols):
                col_wf = [0] * len(MultipleAlignmentArray.alphabet)
                # for the current aligned column, iterate over sequences
                for j in range(self.__nrows):
                    aa = self.msa_arr[j, i]
                    if aa != self.gap:
                        aa_index = MultipleAlignmentArray.alphabet.index(aa)
                        col_wf[aa_index] += sw[j]
                col_wf = [wf + self.__pc * MultipleAlignmentArray.bg_rf[i] for i, wf in enumerate(col_wf)]
                all_wf.append(col_wf)
# the following code is extremely inefficient
#             self.__wf = np.array([[sum((self.msa_arr[:, i] == aa) * self.sequence_weights_henikoff()) 
#                                    for aa in MultipleAlignmentArray.alphabet]
#                                   for i in range(self.__ncols)])
        return np.array(all_wf)
    
    @property
    def icf(self):
        """
            @summary: independent counts
            @see: Pei, 2001, Bioinformatics
        """
        if self.__icf is None:
            self.__icf = []
            for col in range(self.__ncols):
                m = np.array([self.msa_arr[:, col] == aa for aa in MultipleAlignmentArray.alphabet])
                sub_aln = np.array([self.msa_arr[m[i]] for i in range(20) ])
                f = np.array(
                             [sum(
                                  [sum(np.in1d(MultipleAlignmentArray.alphabet, sub_aln[i][:,j])) for j in range(self.__ncols)]                                  
                            ) / float(self.__ncols) for i in range(20)]
                        )
                self.__icf.append(np.log(1.0 - f / 20.0) / np.log(0.95))
        return np.array(self.__icf)  
    
    
    @property
    def gap_f(self):
        """
        """
        if self.__gap_f is None:
            self.__gap_f = np.array([sum(self.msa_arr[:, i] == self.gap) for i in range(self.__ncols)])
        return self.__gap_f
     
    @property
    def gap_rf(self):
        """
        """
        if self.__gap_rf is None:
            self.__gap_rf = np.true_divide(self.gap_f, self.msa_arr.shape[0])
        return self.__gap_rf
        
    
    @property
    def group_f( self ):
        """
        """
        if self.__group_f is None:
            print( "Counting the occurrence of each amino acid group in each aligned column." )
            self.__group_f = np.array( [[sum( np.in1d( self.msa_arr[:, i], g ) ) for g in MultipleAlignmentArray.groups]
                                              for i in range( self.__ncols )] )
        return self.__group_f
    
    @property
    def group_rf(self):
        """
        """
        if self.__group_rf is None:
            print("Calculating the probability that each amino acid group in each aligned column.")
            self.__group_rf = np.true_divide( self.group_f, self.group_f.sum(axis=1)[:, np.newaxis])
        return self.__group_rf
    
    def shannon_entropy( self ):
        """
            @see: Capra, 2007, Bioinformatics
        """
        se = -np.sum( self.rf * np.where( self.rf != 0, np.log2( self.rf ), 0 ), axis = 1 )
        # penalize the conservation score by multiplying it by the fraction of non-gapped positions in the column
        return (np.log2(20) - se) * (1.0 - self.gap_rf)
    
    def shannon_entropy_sequence_weighted(self):
        """
            @summary: 
            @see: Pei, 2001, Bioinformatics
        """
        sesw = -np.sum(self.wf * np.where(self.wf != 0, np.log2(self.wf), 0), axis = 1)
        return (np.log2(20) - sesw) * (1.0 - self.gap_rf)
        
    def shannon_entropy_property_based(self):
        """
            @see: Predicting functionally important residues from sequence conservation, 2007, Bioinformatics
        """
        sepb = -np.sum(self.group_rf * np.where(self.group_rf != 0, np.log2(self.group_rf), 0), axis = 1)
        # penalize the conservation score by multiplying it by the fraction of non-gapped positions in the column
        return (np.log2(6) - sepb) * (1.0 - self.gap_rf)
    
    def relative_entropy(self):
        """
            @summary: the higher the relative entropy from background distribution, the more conserved the aligned position
        """
        re = np.sum(self.rf * np.where(self.rf != 0, np.log2(self.rf / MultipleAlignmentArray.bg_rf), 0), axis = 1)
        return re * (1.0 - self.gap_rf)
    
    from substitution_matrices import BLOSUM50_QIJ
    def von_neumann_entropy(self, rfreq=rf, s_matrix=BLOSUM50_QIJ):
        """
            @see: Predicting functionally important residues from sequence conservation, 2007, Bioinformatics
            and Are protein-protein interfaces more conserved in sequence than the rest of the protein surface?, 2004, Protein Science
        """
        vne = []
        for i in range(self.__ncols):
            w = np.diag(rfreq[i].dot(s_matrix))
            eigv_w = np.linalg.eigvals(w)
            vne.append(-np.sum(eigv_w * np.where(eigv_w != 0, np.log2(eigv_w), 0)))
        return (np.log2(20) - vne) * (1.0 - self.gap_rf)
    
    def variance_based_conservation(self):
        """
            @summary: Conservation reaches its maximum for the position occupied by an
            invariant amino acid whose frequency in the whole alignment is minimal.
            @see: Pei, 2001, Bioinformatics
        """
        # deviation of column frequencies of amino acids from their overall frequencies
        dev = self.rf - np.true_divide(self.f.sum(axis=0), sum(self.f.sum(axis=0)))
        vc = np.sqrt(np.sum( dev**2, axis = 1))
        return vc

    @property
    def aa_in_col( self ):
        """
            @summary: NumPy array, each row of this array indicates which amino acid types are 
            represented in the corresponding aligned column
        """
        if self.__aa_in_col is None:
            print( "Tallying which amino acid types are represented in each aligned column." )
            self.__aa_in_col = np.array( [[aa in self.msa_arr[:, i] for aa in MultipleAlignmentArray.alphabet]
                                         for i in range( self.__ncols )] )
        return self.__aa_in_col
    
    def wu_kabat( self ):
        """
        """
        n = self.f.max( axis = 1 ) 
        aas_in_cols = self.aa_in_col
        k = aas_in_cols.sum( axis = 1 )
        N = self.msa_arr.shape[0]
        return np.true_divide( k, n ) * N
    
    def sequence_weights_henikoff( self ):
        """
            @summary: Each different residue type is awarded an equal share of the weight,
            this equal share of weight is divided equally among the sequences sharing the same
            residue. For each sequence, the contributions from each position are summed to
            give a sequence weight.
            @see: Position-based Sequence Weights, Henikoff, JMB, 1994
        """
        sequence_weights_henikoff = []
        for seq_index in range( self.msa_arr.shape[0] ):
            # ndarray that represents the current aligned sequence
            mask = np.array( MultipleAlignmentArray.alphabet ) == self.msa_arr[seq_index][:, np.newaxis]
            # extract the frequency of each amino acid of the current sequence in each aligned column
            seq_aa_freqs = self.f[mask]
            # some aligned sequences may have gaps, these gapped columns will not be considered
            num_aa_types_in_col = self.aa_in_col.sum( axis = 1 )[mask.sum( axis = 1 ) == 1]
            # add the sequence weight of current sequence
            sequence_weights_henikoff.append( np.sum( 1.0 / ( num_aa_types_in_col * seq_aa_freqs ) ) )
        return np.array( sequence_weights_henikoff ) / sum( sequence_weights_henikoff )

    def pssm( self ):
        """
            @summary: A position-specific scoring matrix M(p,a) is composed of N rows and 21 columns,
            where N is the length of the subject sequence. The row p corresponds to a sequence position
            of the subject. The first 20 columns of each row specify the score for finding, at that position
            in the subject, each of the 20 amino acid residues. An additional column contains a penalty for
            insertions or deletions at that position.
            @details: 
            @see: Profile analysis: Detection of distantly related proteins, 1987, PNAS
        """
        return np.log(self.rf / self.bg_rf)
        
