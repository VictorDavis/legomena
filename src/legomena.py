
# bloody dependencies
from collections import Counter, namedtuple
from scipy.misc import comb as nCr
import numpy as np
import os
import pandas as pd
from sklearn.linear_model import LinearRegression
import zipfile

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
        print("Frequency distribution accessible as <corpus>.fdist")

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
        fdist = Counter(tokens)
        fdist = dict(fdist)

        # return
        return fdist

    #
    def getLegoVectorK(self, fdist:dict = None, relative_to:list = None) -> np.ndarray:
        '''
        Computes the frequency distribution *of* a frequency distribution.
        Ex: k(banana) -> k({a:3, n:2, b:1}) -> {0:0, 1:1, 2:1, 3:1} -> (0, 1, 1, 1)
            - relative_to: Pass "mother" list of tokens for a count for k_0 = number of types *not* in found
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

        # impute k[0] = the number of types *not* found in the distribution
        k[0] = self.N - sum(k)
        # assert len( set(self.tokens) - set(fdist.keys()) ) == k[0] # definition

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
        return [ type for type, freq in self.fdist.items() if freq == n ]

    #
    def hapaxes(self):
        '''Returns list of hapaxes (words that only appear once)'''

        # hapax = 1-legomena
        return self.nlegomena(1)

    #
    def buildTTRCurve(self, seed:int = None, legomena_upto:int = None) -> pd.DataFrame:
        '''
        Samples the corpus at intervals 1/res, 2/res, ... 1 to build a Type-Token Relation
        curve consisting of points (m, n)
            - seed: (optional) Set a random seed for reproduceability
            - legomena_upto: (optional) include n-legomena counts up to given limit
        '''

        # set random seed
        np.random.seed(seed)

        # sample the corpus
        m_choices = [ int((x+1) * self.M / self.res) for x in range(self.res) ]
        TTR = []
        for m_tokens in m_choices:

            # count types & tokens in random sample of size <m_tokens>
            tokens = np.random.choice(self.tokens, m_tokens, replace = False)
            n_types = self.countTypes(tokens)
            row_tuple = [m_tokens, n_types]

            # calculate legomena k-vector expression of tokens sample
            if legomena_upto is not None:
                fdist   = self.getFrequencyDistribution(tokens)
                lego_k  = self.getLegoVectorK(fdist)
                lego_k  = lego_k[:legomena_upto] # only care about n < upto
                pad_width = legomena_upto - len(lego_k)
                if pad_width > 0: # i.e. len(k) < upto
                    lego_k = np.pad(lego_k, (0, pad_width), mode = 'constant') # pad with zeros
                row_tuple += list(lego_k)

            # append row
            row_tuple = tuple(row_tuple)
            TTR.append(row_tuple)

        # save to self.TTR as dataframe
        colnames = ["m_tokens", "n_types"]
        if legomena_upto is not None:
            colnames += [ "lego_" + str(x) for x in range(legomena_upto) ]
        self.TTR = pd.DataFrame(TTR, columns = colnames)

        # return
        print("Type-Token Relation data accessible as <corpus>.TTR")
        return self.TTR

    #
    def fitHeaps(self) -> LinearRegression:
        '''Finds best-fit (K,B) parameters for Heaps Law'''

        # use linear regression to fit K, B to TTR data on a log/log scale
        df = self.TTR
        model = LinearRegression()
        log_m = np.log(df.m_tokens.values)
        log_n = np.log(df.n_types.values)
        X_values = log_m.reshape(-1, 1)
        model.fit(X = X_values, y = log_n)
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

    #
    def transformationMatrix(self, x:float, D:int, mask = None) -> np.ndarray:
        '''Creates a transformation matrix A (of dimension DxD) which operates on k to simulate sampling.'''

        # max dimension 256
        D = min(D, 256)

        # matrix definition
        x = float(x)
        x_ = 1. - x
        vec_range  = np.arange(D)
        mat_horz   = np.meshgrid(vec_range, np.ones(D))[0]
        mat_vert   = np.meshgrid(np.ones(D), np.arange(D))[1]
        mat_nCr    = nCr(mat_horz, mat_vert) # nCr(i, j)
        vec_xrange = x**np.arange(D)
        mat_xvert  = np.meshgrid(np.ones(D), vec_xrange)[1] # x**j
        mat_triu   = np.triu(mat_horz - mat_vert)
        mat_xtriu  = np.triu(x_**mat_triu) # x_**(i-j)
        mat_nCr    = np.nan_to_num(mat_nCr)
        mat_xvert  = np.nan_to_num(mat_xvert)
        mat_xtriu  = np.nan_to_num(mat_xtriu)
        mat_x_x    = mat_xvert * mat_xtriu
        mat_x_x    = np.nan_to_num(mat_x_x)
        A_x        = mat_nCr * mat_x_x
        # A_x = np.array([ nCr(i, j) * x**j * x_**(i-j) for j in range(D) ])
        # A_x = np.nan_to_num(A_x) # ignore nans

        # return A_x matrix
        return A_x

    #
    def transform(self, x:float) -> np.ndarray:
        '''
        Transforms legomena vector k into k' by operating on it by transformation matrix Ax,
        where A(x) = [a_xij] for x the sampling proportion. In English, sampling proportion x
        of a corpus expressible as k results in a smaller corpus expressible as k'
        '''

        # calculate A_x
        D = len(self.k)
        A_x = self.transformationMatrix(x, D)
        D = A_x.shape[0]

        # transform k
        k  = self.k[:D]
        k_ = np.dot(A_x, k)

        # NOTE: max A_x dimension = 256x256 for computational reasons, therefore impute tail with
        #       smooth uniform remainder of expectations
        pad = len(self.k) - D
        err = self.N - sum(k_)
        eps = err / pad
        padding = np.array([ eps ] * pad)
        k_ = np.concatenate((k_, padding))
        assert int(sum(k_)) == self.N

        # return transformed vector
        return k_

    # apply new log formula y = f(x) for proportions y,x in [0,1]
    def log_formula(self, x:float):
        '''
        Applies the normalized sublinear growth x -> y formula [0,1] -> [0,1]
        Types = N_z * log_formula(Tokens / M_z)
        Returns:
            - E_x: Expected proportion of types for given proportion of tokens x
            - k_frac: Expected proportions of n-legomena for n = 0 to 5
        '''

        EPS = 1e-3
        D   = 6
        x   = x.reshape(-1)
        k_frac = np.zeros((len(x), D))
        logx  = np.log(x)
        x2, x3, x4, x5, x6 = [ x**i for i in range(2, D+1) ]
        E_x   = logx*x/(x-1)
        k_frac[:,0] = (x - logx*x - 1)/(x-1)
        k_frac[:,1] = (x2 - logx*x - x)/(x-1)**2
        k_frac[:,2] = (x3 - 2*logx*x2 - x)/2/(x-1)**3
        k_frac[:,3] = (2*x4 - 6*logx*x3 + 3*x3 - 6*x2 + x)/6/(x-1)**4
        k_frac[:,4] = (3*x5 - 12*logx*x4 + 10*x4 - 18*x3 + 6*x2 - x)/12/(x-1)**5
        k_frac[:,5] = (12*x6 - 60*logx*x5 + 65*x5 - 120*x4 + 60*x3 - 20*x2 + 3*x)/60/(x-1)**6

        # return proportions
        return E_x, k_frac

    # predicts number of n-legomena based on input (m) and free parameters (Mz,Nz)
    def predict(self, m_tokens:int, M_z:int = None, N_z:int = None):
        '''
        Applies the log formula to predict number of types as a function of tokens.
            - m_tokens: Number of tokens
            - M_z, N_z: Number of tokens, types of an optimum sample (model parameters)
        Returns:
            - E_m: Expected number of types given m tokens
            - k: (array) Expected n-legomena counts for n = 0 to 5
        '''

        # default to fitted params
        if M_z is None:
            M_z = self.M_z
        if N_z is None:
            N_z = self.N_z

        # scale down to normalized log formula
        x = np.array(m_tokens) / M_z
        E_x, k_frac = self.log_formula(x)
        E_m = np.round(N_z * E_x)
        k   = np.round(N_z * k_frac)

        # de-vectorize if only one point was passed
        if len(E_m) == 1:
            E_m = E_m[0]
            k   = k.reshape(-1)

        # return scaled up predictions
        return E_m, k

    # At what sample size do the n-legomena counts most closely resemble a perfect Zipf Distribution?
    # When 1/ln(x)+1/(1-x) = h_obs for h_obs the observed proportion of hapaxes
    # Note: Because of the nature of this function, Newton's Method is not ideal.
    #       Instead, we use a binary search to find f(x)-h_obs = 0, which
    #       works nicely since f(x) decreases monotonically from 0 to inf
    def inverse_h(self, h_obs):
        '''Solves 1/ln(x)+1/(1-x) = h_obs for x given h_obs'''

        x_prev = 0
        x      = .99
        dx     = .5
        timer  = 1
        EPS    = 1e-8 # ~27 iterations
        while (dx > EPS) & (timer < 999):
            if (x == 0) | (x == 1): # f(x) can't be evaluated @0,1
                x += EPS # jigger x slightly
            fx = 1/np.log(x) + 1/(1-x) # f(x)
            if fx > h_obs:# if f(x) > h_obs then x is too low
                if x + dx == x_prev:
                    dx = dx/2
                x_prev = x
                x += dx
            elif fx < h_obs: # if f(x) < h_obs then x is too high
                if (x - dx == x_prev):
                    dx = dx/2
                while x - dx <= 0: # do not let x go negative
                    dx = dx/2
                x_prev = x
                x -= dx
            else: # if f(x) = h_obs then x is just right
                # print(f"Found x in {timer} iterations")
                return x
            timer += 1
            # print(f"(x, f(x)) = ({x}, {fx})")
        # print(f"Found x = {x} in {timer} iterations.")
        return x

    # find M_z, N_z
    def fit(self):
        '''Uses observed whole-corpus hapax:type ratio to infer M_z, N_z optimum sample size.'''

        # infer z from h_obs
        h_obs = self.k[1] / self.N      # observed hapax frac
        z     = self.inverse_h(h_obs)   # inferred corpus:optimum ratio
        M_z   = self.M / z
        N_z   = self.N * (z-1) / z / np.log(z)
        M_z, N_z = int(M_z), int(N_z)

        # return optimum sample size M_z, N_z
        self.M_z, self.N_z = M_z, N_z
        print("Log model accessible as <corpus>.M_z, .N_z")
        return M_z, N_z

# wrapper class for retrieving corpus from Standard Project Gutenberg Corpus
class SPGC:

    # get corpus by Project Gutenberg book ID
    def get(pgid):

        # extract contents of "counts" text file
        SPGC  = "SPGC-counts-2018-07-18"
        fname = f"{SPGC}/PG{pgid}_counts.txt"
        z = zipfile.ZipFile(f"../../data/{SPGC}.zip")
        z.extract(fname, "./")
        z.close()
        df = pd.read_csv(fname, header = -1, delimiter = "\t")
        os.remove(fname)
        os.rmdir(SPGC)

        # reconstruct "bag of words" from counts
        df.columns = ["word", "freq"]
        nested_list = [ [row.word] * row.freq for idx, row in df.iterrows() ]
        words = [ word for sub_list in nested_list for word in sub_list ]
        corpus = Corpus(words)

        # return corpus object
        return corpus
