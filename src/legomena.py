
# bloody dependencies
import collections
import numpy as np
import pandas as pd

# main class
class Legomena:

    ## object properties
    tokens = None           # corpus under study, list of strings
    fdist  = None           # frequency distribution of tokens in corpus
    k      = None           # k-vector of n-legomena counts
    M      = None           # number of tokens in the corpus
    N      = None           # number of types in the corpus
    res    = 100            # number of samples to take when building a TTR curve
    ttr_df = None           # dataframe containing type/token counts from corpus samples

    #
    def __init__(self, tokens:list):
        '''Initializes handler class with a corpus as a list of tokens.'''

        # store passed tokens & pre-calculate a few items
        self.tokens = tokens
        self.M = self.countTokens()
        self.N = self.countTypes()
        self.fdist = self.getFrequencyDistribution()
        self.k = self.getLegoVectorK()

        # report
        print("Number of tokens (self.M):", self.M)
        print("Number of types  (self.N):", self.N)
        print("Legomena vector  (self.k):", self.k)

    #
    def getFrequencyDistribution(self, tokens:list = None, normalized:bool = False) -> dict:
        '''Takes a list of tokens and returns a dictionary of {token: count} pairs.'''

        # pre-calculated frequency distribution
        if tokens is None and self.fdist is not None:
            return self.fdist

        # pre-calculate based on initially defined corpus
        if tokens is None:
            tokens = self.tokens

        # calculate frequency distribution
        fdist = collections.Counter(self.tokens)
        fdist = dict(fdist)
        return fdist

    #
    def getLegoVectorK(self, fdist:dict = None) -> np.ndarray:
        '''
        Computes the frequency distribution *of* a frequency distribution.
        Ex: k(banana) -> {a:3, n:2, b:1} -> {0:0, 1:1, 2:1, 3:1} -> (0, 1, 1, 1)
        '''

        # pre-calculated k-vector of corpus
        if fdist is None and self.k is not None:
            return self.k

        # pre-calculate based on previously calculated fdist
        if fdist is None:
            fdist = self.fdist

        # encodes frequency distribution as a vector of n-legomena counts
        fdist_vals = [ val for key, val in fdist.items() ]
        fdist_n, fdist_kn = np.unique(fdist_vals, return_counts = True)

        # impute zeros for elements of k for which no words occur n times, including n=0
        df = pd.DataFrame({"n": fdist_n, "kn": fdist_kn})
        max_n = df.n.max()
        df.set_index("n", inplace = True)
        df_complete = pd.DataFrame(index = range(max_n + 1))
        df = pd.concat([df_complete, df], axis = 1)
        df.fillna(0, inplace = True)
        k = df.astype({"kn":int}).kn.values

        # return vector
        return k

    #
    def countTokens(self, tokens:list = None):
        '''Returns the number of tokens in the corpus.'''

        # pre-calculated M
        if tokens is None and self.M is not None:
            return self.M

        # calculate M of self-defined corpus
        if tokens is None:
            tokens = self.tokens

        # size of corpus
        return len(tokens)

    #
    def countTypes(self, tokens:list = None):
        '''Returns the number of types in the corpus.'''

        # pre-calculated N
        if tokens is None and self.N is not None:
            return self.N

        # calculate N of self-defined corpus
        if tokens is None:
            tokens = self.tokens

        # size of set of unique elements
        # NOTE: np.unique() faster than set(x)
        return len(np.unique(tokens))

    #
    def minicorpus(self, n:int):
        '''Returns (kn,k)-minicorpus, list of all words occurring exactly n times.'''

        # list of n-legomena
        nlegomena = self.nlegomena(n)

        # all occurrences of n-legomena
        return [ token for token in tokens if token in nlegomena ]

    #
    def nlegomena(self, n:int):
        '''Returns list of tokens occurring exactly n times in the corpus.'''

        # all times occurring exactly n times
        return [ type for type, frequency in self.fdist.items() if freq == 1 ]

    #
    def fitHeaps(self) -> tuple:
        '''Finds best-fit K,B parameters for Heaps Law'''
        return None
