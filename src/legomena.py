
# bloody dependencies
from collections import Counter, namedtuple
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# main class
class Corpus:

    ## object properties
    tokens = None           # corpus under study, list of strings
    fdist  = None           # frequency distribution of tokens in corpus
    k      = None           # k-vector of n-legomena counts
    M      = None           # number of tokens in the corpus
    N      = None           # number of types in the corpus
    res    = 100            # number of samples to take when building a TTR curve
    TTR    = None           # dataframe containing type/token counts from corpus samples
    heaps_K = None          # best-fit Heap's coefficient N = KM^B
    heaps_B = None          # best-fit Heap's exponent N = KM^B

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
        print("Number of tokens (<corpus>.M):", self.M)
        print("Number of types  (<corpus>.N):", self.N)
        print("Legomena vector  (<corpus>.k):", self.k[:9])

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
        fdist = Counter(self.tokens)
        fdist = dict(fdist)

        # return
        print("Frequency distribution accessible as <corpus>.fdist")
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
        k = np.array(k)

        # return vector
        print("Legomena counts accessible as <corpus>.k")
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
        return [ type for type, freq in self.fdist.items() if freq == n ]

    #
    def hapaxes(self):
        '''Returns list of hapaxes (words that only appear once)'''

        # hapax = 1-legomena
        return self.nlegomena(1)

    #
    def buildTTRCurve(self, seed:int = None) -> pd.DataFrame:
        '''
        Samples the corpus at intervals 1/res, 2/res, ... 1 to build a Type-Token Relation
        curve consisting of points (m, n)
            - seed: (optional) Set a random seed for reproduceability
        '''

        # set random seed
        np.random.seed(seed)

        # sample the corpus
        m_choices = [ int((x+1) * self.M / self.res) for x in range(self.res) ]
        TTR = []
        for m_tokens in m_choices:
            tokens = np.random.choice(self.tokens, m_tokens, replace = False)
            n_types = self.countTypes(tokens)
            TTR.append((m_tokens, n_types))

        # save to self.TTR as dataframe
        self.TTR = pd.DataFrame(TTR, columns = ["m_tokens", "n_types"])

        # annotate for log-log plotting & for fitting Heap's Law
        self.TTR["log_m"] = np.log(self.TTR.m_tokens)
        self.TTR["log_n"] = np.log(self.TTR.n_types)

        # return
        print("Type-Token Relation curve accessible as <corpus>.TTR")
        return self.TTR

    #
    def fitHeaps(self) -> LinearRegression:
        '''Finds best-fit (K,B) parameters for Heaps Law'''

        # use linear regression to fit K, B to TTR data on a log/log scale
        df = self.TTR
        model = LinearRegression()
        X_values = df.log_m.values.reshape(-1, 1)
        model.fit(X = X_values, y = df.log_n)
        B = model.coef_[0]
        K = np.exp(model.intercept_)
        self.heaps_K, self.heaps_B = K, B

        # return parameters
        print("Best-fit Heap's Law parameters E(m) = Km^B (K,B) = ", (K, B))
        print("Heap's Law model accessible as <corpus>.heaps(m)")
        return (K, B)

    #
    def heaps(self, m:int) -> int:
        '''Runs best-fit Heap's Law model to predict n_types = E(m_tokens)'''

        # assumes model has been fitted
        K, B = self.heaps_K, self.heaps_B
        E_m = np.round(K * np.power(m, B))

        # return
        return E_m

    # infinite series model with coefficients determined by legomena vector k
    # eqn (8b) in https://arxiv.org/pdf/1901.00521.pdf
    def iseries(self, m:int) -> int:
        '''Runs infinite series model to predict N = E(M)'''

        # sum over legomena coefficients k
        m = np.array(m)
        x = m / self.M
        exponents = range(len(self.k))
        terms = np.array([ np.power(1 - x, n) for n in exponents ])
        k_0 = np.dot(self.k, terms)
        E_x = np.round(self.N - k_0)

        # return
        return E_x
