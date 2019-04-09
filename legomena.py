# bloody dependencies
from collections import Counter, namedtuple
from scipy.special import comb as nCr
import numpy as np
import os
import pandas as pd
from scipy.optimize import fsolve, curve_fit
from sklearn.linear_model import LinearRegression
import zipfile

# environment vars
DATAPATH = os.getenv("DATAPATH", "data")

# main class
class Corpus:

    ## object properties
    tokens = None  # corpus under study, list of strings
    types = None  # lexicon, ranked order
    fdist = None  # frequency distribution of tokens in corpus
    k = None  # k-vector of n-legomena counts
    M = None  # number of tokens in the corpus
    N = None  # number of types in the corpus
    res = 100  # number of samples to take when building a TTR curve
    TTR = None  # dataframe containing type/token counts from corpus samples
    heaps_K = None  # best-fit Heap's coefficient N = KM^B
    heaps_B = None  # best-fit Heap's exponent N = KM^B
    _mat = None  # internal cache for speeding up .transform() function

    #
    def __init__(self, tokens: list, verbose: bool = True):
        """Initializes handler class with a corpus as a list of tokens."""

        # store passed tokens & pre-calculate a few items
        self.tokens = tokens
        self.M = self.countTokens()
        self.N = self.countTypes()
        self.fdist = self.getFrequencyDistribution()
        self.types = self.fdist.type.values
        self.k = self.getLegoVectorK()

        # report
        if verbose:
            print("Number of tokens (<corpus>.M):", self.M)
            print("Number of types  (<corpus>.N):", self.N)
            print("Legomena vector  (<corpus>.k):", self.k[:9])
            print("Frequency distribution accessible as <corpus>.fdist")

    #
    def getFrequencyDistribution(
        self, tokens: list = None, normalized: bool = False
    ) -> pd.core.frame.DataFrame:
        """Takes a list of tokens and returns a dataframe of type/rank/frequency counts."""

        # pre-calculated frequency distribution
        if tokens is None and self.fdist is not None:
            return self.fdist

        # pre-calculate based on initially defined corpus
        if tokens is None:
            tokens = self.tokens

        # calculate frequency distribution
        fdist = Counter(tokens)
        fdist = dict(fdist)
        fdist = pd.DataFrame.from_dict(fdist, orient="index").reset_index()
        fdist.columns = ["type", "freq"]
        fdist.sort_values("freq", ascending=False, inplace=True)
        fdist["rank"] = range(1, fdist.shape[0] + 1)

        # return
        return fdist

    #
    def getLegoVectorK(
        self, fdist: pd.core.frame.DataFrame = None, relative_to: list = None
    ) -> np.ndarray:
        """
        Computes the frequency distribution *of* a frequency distribution.
        Ex: k(banana) -> k({a:3, n:2, b:1}) -> {0:0, 1:1, 2:1, 3:1} -> (0, 1, 1, 1)
            - relative_to: Pass "mother" list of tokens for a count for k_0 = number of types *not* in found
        """

        # pre-calculated k-vector of corpus
        if fdist is None and self.k is not None:
            return self.k

        # pre-calculate based on previously calculated fdist
        if fdist is None:
            fdist = self.fdist

        # impute zeros for elements of k for which no words occur n times, including n=0
        df = fdist.freq.value_counts().to_frame()
        max_freq = df.index.max()
        df_complete = pd.DataFrame(index=range(max_freq + 1))
        df = pd.concat([df_complete, df], axis=1)
        df.fillna(0, inplace=True)
        k = df.astype({"freq": int}).freq.values
        k = np.array(k)

        # impute k[0] = the number of types *not* found in the distribution
        k[0] = self.N - sum(k)

        # return vector
        return k

    #
    def countTokens(self, tokens: list = None):
        """Returns the number of tokens in the corpus."""

        # pre-calculated M
        if tokens is None and self.M is not None:
            return self.M

        # calculate M of self-defined corpus
        if tokens is None:
            tokens = self.tokens

        # size of corpus
        return len(tokens)

    #
    def countTypes(self, tokens: list = None):
        """Returns the number of types in the corpus."""

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
    def minicorpus(self, n: int):
        """Returns (kn,k)-minicorpus, list of all words occurring exactly n times."""

        # list of n-legomena
        nlegomena = self.nlegomena(n)

        # all occurrences of n-legomena
        return [token for token in tokens if token in nlegomena]

    #
    def nlegomena(self, n: int):
        """Returns list of tokens occurring exactly n times in the corpus."""

        # all times occurring exactly n times
        df = self.fdist
        return df[df.freq == n].type.values

    #
    def hapaxes(self):
        """Returns list of hapaxes (words that only appear once)"""

        # hapax = 1-legomena
        return self.nlegomena(1)

    #
    def buildTTRCurve(
        self, seed: int = None, legomena_upto: int = None
    ) -> pd.DataFrame:
        """
        Samples the corpus at intervals 1/res, 2/res, ... 1 to build a Type-Token Relation
        curve consisting of points (m, n)
            - seed: (optional) Set a random seed for reproduceability
            - legomena_upto: (optional) include n-legomena counts up to given limit
        """

        # set random seed
        np.random.seed(seed)

        # sample the corpus
        m_choices = [int((x + 1) * self.M / self.res) for x in range(self.res)]
        TTR = []
        for m_tokens in m_choices:

            # count types & tokens in random sample of size <m_tokens>
            tokens = np.random.choice(self.tokens, m_tokens, replace=False)
            n_types = self.countTypes(tokens)
            row_tuple = [m_tokens, n_types]

            # calculate legomena k-vector expression of tokens sample
            if legomena_upto is not None:
                fdist = self.getFrequencyDistribution(tokens)
                lego_k = self.getLegoVectorK(fdist)
                rank1 = len(lego_k)  # frequency of most frequent token
                lego_k = lego_k[:legomena_upto]  # only care about n < upto
                pad_width = legomena_upto - len(lego_k)
                if pad_width > 0:  # i.e. len(k) < upto
                    lego_k = np.pad(
                        lego_k, (0, pad_width), mode="constant"
                    )  # pad with zeros
                row_tuple += list(lego_k) + [rank1]

            # append row
            row_tuple = tuple(row_tuple)
            TTR.append(row_tuple)

        # save to self.TTR as dataframe
        colnames = ["m_tokens", "n_types"]
        if legomena_upto is not None:
            colnames += ["lego_" + str(x) for x in range(legomena_upto)]
            colnames += ["rank1"]
        self.TTR = pd.DataFrame(TTR, columns=colnames)

        # return
        print("Type-Token Relation data accessible as <corpus>.TTR")
        return self.TTR

    #
    def fitHeaps(self) -> LinearRegression:
        """Finds best-fit (K,B) parameters for Heaps Law"""

        # use linear regression to fit K, B to TTR data on a log/log scale
        df = self.TTR
        model = LinearRegression()
        log_m = np.log(df.m_tokens.values)
        log_n = np.log(df.n_types.values)
        X_values = log_m.reshape(-1, 1)
        model.fit(X=X_values, y=log_n)
        B = model.coef_[0]
        K = np.exp(model.intercept_)
        self.heaps_K, self.heaps_B = K, B

        # return parameters
        print("Best-fit Heap's Law parameters E(m) = Km^B (K,B) = ", (K, B))
        print("Heap's Law model accessible as <corpus>.heaps(m)")
        return (K, B)

    #
    def heaps(self, m: int) -> int:
        """Runs best-fit Heap's Law model to predict n_types = E(m_tokens)"""

        # assumes model has been fitted
        K, B = self.heaps_K, self.heaps_B
        E_m = np.round(K * np.power(m, B))

        # return
        return E_m

    # infinite series model with coefficients determined by legomena vector k
    # eqn (8b) in https://arxiv.org/pdf/1901.00521.pdf
    def iseries(self, m: int) -> int:
        """Runs infinite series model to predict N = E(M)"""

        # sum over legomena coefficients k
        m = np.array(m)
        x = m / self.M
        exponents = range(len(self.k))
        terms = np.array([np.power(1 - x, n) for n in exponents])
        k_0 = np.dot(self.k, terms)
        E_x = np.round(self.N - k_0)

        # return
        return E_x

    #
    def transformationMatrix(self, x: float, D: int, mask=None) -> np.ndarray:
        """
        Creates a transformation matrix A (of dimension DxD) which operates on k to simulate sampling.
            x: float between 0 and 1
            D: desired dimension of matrix A = DxD
            mask: (optional) to save memory (~10x, depending), ignore matrix elements corresponding to zeros in the k-vector
        """

        # cache elements of the calculation invariant to x
        if (self._mat is None) or (D != len(self.k)):

            # horizontal/vertical indices
            mat_horz = np.repeat(np.arange(D).reshape(1, D), D, axis=0)  # i
            mat_vert = np.repeat(np.arange(D).reshape(D, 1), D, axis=1)  # j
            mat_triu = mat_horz - mat_vert  # i-j

            # coerce into upper triangular
            mat_vert = np.triu(mat_vert)
            mat_triu = np.triu(mat_triu)

            # now reduce matrix considerably along input axis for zero-elements of k-vector
            if mask is not None:
                mat_horz = mat_horz[:, mask]
                mat_vert = mat_vert[:, mask]
                mat_triu = mat_triu[:, mask]

            # compute pascal's matrix
            mat_pasc = nCr(mat_horz, mat_vert)  # nCr(i, j)
            if mask is None:
                mat_pasc = np.triu(mat_pasc)
            else:
                mat_pasc_sq = np.zeros((D, D))
                mat_pasc_sq[:, mask] = mat_pasc
                mat_pasc_sq = np.triu(mat_pasc_sq)
                mat_pasc = mat_pasc_sq[:, mask]

            # assume infinite values will be zeroed out -- TODO: check?
            mat_pasc[np.isinf(mat_pasc)] = 0
            assert np.isinf(mat_pasc).sum() == 0

            # cache
            if D == len(self.k):
                self._mat = {"vert": mat_vert, "triu": mat_triu, "pasc": mat_pasc}
        else:
            mat_vert = self._mat["vert"]
            mat_triu = self._mat["triu"]
            mat_pasc = self._mat["pasc"]

        # verify dimensions
        D_ = D if mask is None else sum(mask)
        assert mat_vert.shape == (D, D_)
        assert mat_triu.shape == (D, D_)
        assert mat_pasc.shape == (D, D_)

        # add x into the recipe
        x = np.array(x)
        x_ = 1.0 - x

        # support 1-dim argument for x
        xdim = np.ndim(x)
        assert xdim in [0, 1]  # scalar or 1-D array
        if xdim == 1:
            R = x.shape[0]
            mat_vert = np.repeat(mat_vert[:, :, np.newaxis], R, axis=2)
            mat_triu = np.repeat(mat_triu[:, :, np.newaxis], R, axis=2)
            mat_pasc = np.repeat(mat_pasc[:, :, np.newaxis], R, axis=2)

        # calculate exponents and combine
        mat_xvert = np.power(x, mat_vert)  # x^j
        mat_xtriu = np.power(x_, mat_triu)  # x_^(i-j)
        mat_xprod = mat_xvert * mat_xtriu  # x^j * x_^(i-j)
        A_x = mat_pasc * mat_xprod  # nCr(i, j) * x^j * x_^(i-j)

        # coerce matrix shape to ( # of x elements, # of k elements, # of non-zero k elements )
        if xdim == 1:
            A_x = np.moveaxis(A_x, 2, 0)
            assert A_x.shape == (R, D, D_)
        else:
            assert A_x.shape == (D, D_)

        # return
        return A_x

    #
    def transform(self, x: float) -> np.ndarray:
        """
        Transforms legomena vector k into k' by operating on it by transformation matrix Ax,
        where A(x) = [a_xij] for x the sampling proportion. In English, sampling proportion x
        of a corpus expressible as k results in a smaller corpus expressible as k'
        """

        # calculate A_x
        D = len(self.k)
        mask = self.k > 0
        A_x = self.transformationMatrix(x, D, mask)

        # transform k
        k = self.k[mask]
        k_ = np.dot(A_x, k)

        # return transformed vector
        return k_

    # apply new log formula y = f(x) for proportions y,x in [0,1]
    def log_formula(self, x: float):
        """
        Applies the normalized sublinear growth x -> y formula [0,1] -> [0,1]
        Types = N_z * log_formula(Tokens / M_z)
        Returns:
            - E_x: Expected proportion of types for given proportion of tokens x
            - k_frac: Expected proportions of n-legomena for n = 0 to 5
        """

        EPS = 1e-3
        D = 6
        x = np.array(x).reshape(-1)

        # initialize return vectors
        _E_x = np.zeros((len(x)))
        _k_frac = np.zeros((len(x), D))

        # two singularities @ x=0 & x=1
        x_iszero = np.abs(x) < EPS  # singularity @ x=0
        x_isone = np.abs(x - 1.0) < EPS  # singularity @ x=1
        x_isokay = ~np.logical_or(x_iszero, x_isone)
        x = x[x_isokay]

        # initialize results
        k_frac = np.zeros((len(x), D))
        logx = np.log(x)
        x2, x3, x4, x5, x6 = [x ** i for i in range(2, D + 1)]
        E_x = logx * x / (x - 1)
        k_frac[:, 0] = (x - logx * x - 1) / (x - 1)
        k_frac[:, 1] = (x2 - logx * x - x) / (x - 1) ** 2
        k_frac[:, 2] = (x3 - 2 * logx * x2 - x) / 2 / (x - 1) ** 3
        k_frac[:, 3] = (2 * x4 - 6 * logx * x3 + 3 * x3 - 6 * x2 + x) / 6 / (x - 1) ** 4
        k_frac[:, 4] = (
            (3 * x5 - 12 * logx * x4 + 10 * x4 - 18 * x3 + 6 * x2 - x)
            / 12
            / (x - 1) ** 5
        )
        k_frac[:, 5] = (
            (12 * x6 - 60 * logx * x5 + 65 * x5 - 120 * x4 + 60 * x3 - 20 * x2 + 3 * x)
            / 60
            / (x - 1) ** 6
        )

        # limiting values @ x=0 & x=1
        _E_x[x_iszero] = 0.0
        _k_frac[x_iszero, :] = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        _E_x[x_isone] = 1.0
        _k_frac[x_isone, :] = np.array(
            [0.0, 1.0 / 2, 1.0 / 6, 1.0 / 12, 1.0 / 20, 1.0 / 30]
        )
        _E_x[x_isokay] = E_x
        _k_frac[x_isokay] = k_frac

        # return proportions
        return _E_x, _k_frac

    # predicts number of n-legomena based on input (m) and free parameters (Mz,Nz)
    def predict(self, m_tokens: int, M_z: int = None, N_z: int = None):
        """
        Applies the log formula to predict number of types as a function of tokens.
            - m_tokens: Number of tokens
            - M_z, N_z: Number of tokens, types of an optimum sample (model parameters)
        Returns:
            - E_m: Expected number of types given m tokens
            - k: (array) Expected n-legomena counts for n = 0 to 5
        """

        # default to fitted params
        if M_z is None:
            M_z = self.M_z
        if N_z is None:
            N_z = self.N_z

        # scale down to normalized log formula
        x = np.array(m_tokens) / M_z
        E_x, k_frac = self.log_formula(x)
        E_m = np.round(N_z * E_x)
        k = np.round(N_z * k_frac)

        # de-vectorize if only one point was passed
        if len(E_m) == 1:
            E_m = E_m[0]
            k = k.reshape(-1)

        # return scaled up predictions
        return E_m, k

    # find M_z, N_z
    def fit(self, optimize: bool = False):
        """
        Uses observed whole-corpus hapax:type ratio to infer M_z, N_z optimum sample size.
            - optimize: (optional) Uses scipy.optimize.curve_fit() to refine initial fit.
        """

        # infer z from h_obs
        h_obs = self.k[1] / self.N  # observed hapax frac
        func = lambda x: 1.0 / np.log(x) + 1.0 / (1.0 - x) - h_obs
        z0 = 0.5
        z = fsolve(func, z0)[0]
        M_z = self.M / z
        N_z = self.N * (z - 1) / z / np.log(z)
        M_z, N_z = int(M_z), int(N_z)
        self.M_z, self.N_z = M_z, N_z

        # minimize MSE on random perturbations of M_z, N_z
        if optimize:
            if self.TTR is None:
                print("To optimize, call <corpus>.buildTTRCurve() first.")
            else:
                func = (
                    lambda m, M_z, N_z: N_z * np.log(m / M_z) * m / M_z / (m / M_z - 1)
                )
                xdata = self.TTR.m_tokens
                ydata = self.TTR.n_types
                params, _ = curve_fit(func, xdata, ydata, p0=(M_z, N_z))
                M_z, N_z = int(params[0]), int(params[1])
                self.M_z, self.N_z = M_z, N_z

        # return optimum sample size M_z, N_z
        print("Log model accessible as <corpus>.M_z, .N_z")
        return M_z, N_z


# wrapper class for retrieving corpus from Standard Project Gutenberg Corpus
class SPGC:
    """
    Standard Project Gutenberg Corpus
    ---------------------------------
    Created by Font-Clos, F. & Gerlach, M. https://arxiv.org/abs/1812.08092
    Snapshot downloaded from: https://zenodo.org/record/2422561/files/SPGC-counts-2018-07-18.zip
    """

    # get corpus by Project Gutenberg book ID
    def get(pgid: int) -> Corpus:
        """
        Retrieves word frequency distribution for book by PGID, reconstructs (unordered) "bag of words" from counts
        Returns Corpus object
        """

        # extract contents of "counts" text file
        SPGC = "SPGC-counts-2018-07-18"
        zname = f"{DATAPATH}/{SPGC}.zip"
        fname = f"{SPGC}/PG{pgid}_counts.txt"
        fobj = None
        try:
            z = zipfile.ZipFile(zname)
            fobj = z.open(fname)
            z.close()
        except Exception as e:
            print(e)  # print but do not raise

        # check for text file
        if fobj is None:
            try:
                fobj = open(f"{DATAPATH}/{fname}")
            except Exception as e:
                print(e)
                raise (e)

        # reconstruct "bag of words" from counts
        df = pd.read_csv(fobj, header=-1, delimiter="\t")
        df.columns = ["word", "freq"]
        nested_list = [[row.word] * row.freq for idx, row in df.iterrows()]
        words = [word for sub_list in nested_list for word in sub_list]
        corpus = Corpus(words)

        # return corpus object
        return corpus

    # get SPGC metadata from csv
    def getMeta():
        """Retrieves metadata.csv and returns as dataframe."""

        # read csv
        df = pd.read_csv(f"{DATAPATH}/SPGC-metadata-2018-07-18.csv")

        # strip out "sound" entries
        df = df[df.type == "Text"]
        assert df.shape == (56466, 9)

        # key records by numeric PG ID
        df.id = df.id.str.replace("PG", "")
        df = df.set_index("id")

        # return
        return df
