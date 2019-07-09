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


class HeapsModel:
    """Heap's Law: Types = K*Tokens^B"""

    # params
    HeapsParams = namedtuple("HeapsParams", ("K", "B"))
    _params = None

    @property
    def params(self) -> HeapsParams:
        assert (
            self._params is not None
        ), "Please use .fit(m_tokens, n_types) to fit the model."
        return self._params

    @params.setter
    def params(self, params_: tuple):
        self._params = self.HeapsParams(*params_)

    @property
    def K(self) -> int:
        """Heap's Law coefficient K."""
        return self.params.K

    @property
    def B(self) -> int:
        """Heap's Law exponent B."""
        return self.params.B

    def fit(self, m_tokens: np.ndarray, n_types: np.ndarray) -> HeapsParams:
        """
        Naive fit: Linear regression on logN = B*logM + expK
        :param m_tokens: Number of tokens, list-like independent variable
        :param n_types: Number of types, list-like dependent variable
        :returns: Heap's Law parameters, as tuple
        """

        # convert to log/log
        log_m = np.log(m_tokens).reshape(-1, 1)
        log_n = np.log(n_types)
        model = LinearRegression()
        model.fit(X=log_m, y=log_n)
        K = np.exp(model.intercept_)
        B = model.coef_[0]
        self.params = (K, B)
        return self.params

    def predict(self, m_tokens: np.ndarray) -> np.ndarray:
        """
        Calculate & return n_types = E(m_tokens) = Km**B
        :param m_tokens: Number of tokens, list-like independent variable
        :returns: Number of types as predicted by Heap's Law
        """

        # allow scalar
        return_scalar = np.isscalar(m_tokens)

        # run model
        K, B = self.params
        n_types = np.round(K * np.power(m_tokens, B))

        # allow scalar
        if return_scalar:
            n_types = n_types.squeeze()

        # return
        return n_types


class KTransformer:
    """Matrix operator which transforms k(x) = A_x*k(1)"""

    @classmethod
    def A_x(cls, x: float, dim: int, mask=None) -> np.ndarray:
        """
        Creates a transformation matrix A (of dimension dim^2) which operates on k to simulate sampling.
        NOTE: Eqn (12) in https://arxiv.org/pdf/1901.00521.pdf
        :param x: float between 0 and 1
        :param dim: desired dimension of matrix A = dim^2
        :param mask: (optional) to save memory (~10x, depending), ignore matrix elements corresponding to zeros in the k-vector
        :returns: Transformation matrix A_x
        """

        # horizontal/vertical indices
        mat_horz = np.repeat(np.arange(dim).reshape(1, dim), dim, axis=0)  # i
        mat_vert = np.repeat(np.arange(dim).reshape(dim, 1), dim, axis=1)  # j
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
            mat_pasc_sq = np.zeros((dim, dim))
            mat_pasc_sq[:, mask] = mat_pasc
            mat_pasc_sq = np.triu(mat_pasc_sq)
            mat_pasc = mat_pasc_sq[:, mask]

        # assume infinite values will be zeroed out -- TODO: check?
        mat_pasc[np.isinf(mat_pasc)] = 0
        assert np.isinf(mat_pasc).sum() == 0

        # verify dimensions
        dim_ = dim if mask is None else sum(mask)
        assert mat_vert.shape == (dim, dim_)
        assert mat_triu.shape == (dim, dim_)
        assert mat_pasc.shape == (dim, dim_)

        # add x into the recipe
        x = np.array(x)
        x_ = 1.0 - x

        # support 1-dim argument for x
        xdim = np.ndim(x)
        assert xdim in [0, 1]  # scalar or 1-dim array
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
            assert A_x.shape == (R, dim, dim_)
        else:
            assert A_x.shape == (dim, dim_)

        # return
        return A_x

    @classmethod
    def transform(cls, k: np.ndarray, x: float) -> np.ndarray:
        """
        Transforms legomena vector k into k' by operating on it by transformation matrix A_x,
        where A_x = [a_xij] for x the sampling proportion. In English, sampling proportion x
        of a corpus expressible as k results in a smaller corpus expressible as k'
        :param k: The k-vector to transform
        :param x: The matrix parameter for A_x
        :returns: The transformed k-vector
        """

        # calculate A_x
        mask = k > 0
        dim = len(k)
        A_x = cls.A_x(x, dim, mask)

        # transform k
        k = k[mask]
        k_ = np.dot(A_x, k)

        # return transformed vector
        return k_


class Corpus(Counter):
    """Wrapper class for logarithmic competitor model to Heap's Law"""

    # object properties
    _k = None  # k-vector of n-legomena counts
    _TTR = None  # dataframe containing type/token counts from corpus samples

    # params
    LogParams = namedtuple("LogParams", ("M_z", "N_z"))
    _params = None

    # user options
    UserOptions = namedtuple(
        "UserOptions", ("resolution", "dimension", "seed", "epsilon")
    )
    _options = UserOptions(
        resolution=100,  # number of samples to take when building a TTR curve
        dimension=7,  # vector size for holding n-legomena counts (zero-based)
        seed=None,  # random number seed for taking samples when building TTR
        epsilon=1e-3,  # proximity to singularity to abort log formula calculation
    )

    @property
    def tokens(self) -> list:
        """The bag of words, list-like of elements of any type."""
        tokens_ = list(self.elements())
        return tokens_

    @property
    def types(self) -> list:
        """The lexicon, ranked by frequency."""
        fdist = self.fdist  # ranked order
        types_ = list(fdist.type.values)
        return types_

    @property
    def fdist(self) -> pd.DataFrame:
        """The frequency distribution, as a dataframe."""
        df = pd.DataFrame.from_dict(self, orient="index").reset_index()
        df.columns = ["type", "freq"]
        df = df.sort_values("freq", ascending=False).reset_index(drop=True)
        return df

    @property
    def k(self) -> np.ndarray:
        """The vector describing the frequency of n-legomena."""
        if self._k is None:
            self._k = Counter(self.values())

        # return as array
        k = self._k
        kmax = max(k.keys())
        karr = np.array([k[i] for i in range(kmax + 1)])
        return karr

    @property
    def M(self) -> int:
        """The number of tokens in the corpus."""
        m_tokens = sum(self.values())
        return m_tokens

    @property
    def N(self) -> int:
        """The number of types in the corpus."""
        n_types = len(self)
        return n_types

    @property
    def options(self) -> UserOptions:
        """Misc low-impact options when computing the TTR curve."""
        return self._options

    @options.setter
    def options(self, opt_: tuple):
        self._options = self.UserOptions(*opt_)
        self._TTR = None
        self._params = None

    @property
    def resolution(self) -> int:
        """The desired resolution (number of samples) of the TTR curve."""
        return self.options.resolution

    @resolution.setter
    def resolution(self, res_: int):
        res, dim, seed, eps = self.options
        res = res_
        self.options = (res, dim, seed, eps)

    @property
    def dimension(self) -> int:
        """The desired dimension of the problem, highest n for which n-legomena counts are included in TTR."""
        return self.options.dimension

    @dimension.setter
    def dimension(self, dim_: int):
        res, dim, seed, eps = self.options
        dim = dim_
        self.options = (res, dim, seed, eps)

    @property
    def seed(self) -> int:
        """Random number seed."""
        return self.options.seed

    @seed.setter
    def seed(self, seed_: int):
        res, dim, seed, eps = self.options
        seed = seed_
        self.options = (res, dim, seed, eps)

    @property
    def epsilon(self) -> float:
        """Distance to singularity to abort log formula calculation."""
        return self.options.epsilon

    @epsilon.setter
    def epsilon(self, eps_: float):
        res, dim, seed, eps = self.options
        eps = eps_
        self.options = (res, dim, seed, eps)

    @property
    def TTR(self) -> pd.DataFrame:
        """DataFrame of type-token relation data."""
        if self._TTR is None:
            self._TTR = self._compute_TTR()
        return self._TTR

    @property
    def params(self) -> LogParams:
        """Fitting parameters M_z, N_z of the logarithmic model."""
        if self._params is None:
            self._params = self.fit()
        return self._params

    @params.setter
    def params(self, params_: tuple):
        self._params = self.LogParams(*params_)

    @property
    def M_z(self) -> int:
        """The number of tokens in this corpus's optimum sample."""
        return self.params.M_z

    @property
    def N_z(self) -> int:
        """The number of types in this corpus's optimum sample."""
        return self.params.N_z

    def summary(self):
        """Print some basic information about this corpus."""
        print("Number of tokens (<corpus>.M):", self.M)
        print("Number of types  (<corpus>.N):", self.N)
        print("Legomena vector  (<corpus>.k):", self.k[:9])

    #
    def nlegomena(self, n: int):
        """List of types occurring exactly n times in the corpus."""
        nlegomena_ = [typ for typ, freq in self.items() if freq == n]
        return nlegomena_

    @property
    def hapax(self):
        """List of hapax (words that appear exactly once)"""
        return self.nlegomena(1)

    @property
    def dis(self):
        """List of dis legomena (words that appear exactly twice)"""
        return self.nlegomena(2)

    @property
    def tris(self):
        """List of tris legomena (words that appear exactly three times)"""
        return self.nlegomena(3)

    @property
    def tetrakis(self):
        """List of tetrakis legomena (words that appear exactly four times)"""
        return self.nlegomena(4)

    @property
    def pentakis(self):
        """List of pentakis legomena (words that appear exactly five times)"""
        return self.nlegomena(5)

    #
    def _compute_TTR(self) -> pd.DataFrame:
        """
        Samples the corpus at intervals 1/res, 2/res, ... 1 to build a Type-Token Relation
        curve consisting of <resolution> points <(m, n)>
        :returns: Dataframe of token/type/legomena counts
        """

        # user options
        res, dim, seed, eps = self.options

        # retrieve/reconstruct tokens
        tokens = self.tokens

        # set random seed
        np.random.seed(seed)

        # sample the corpus
        m_choices = [int((x + 1) * self.M / res) for x in range(res)]
        TTR = []
        for m_tokens in m_choices:

            # count types & tokens in random sample of size <m_tokens>
            # NOTE: WITHOUT REPLACEMENT
            subtokens = np.random.choice(tokens, m_tokens, replace=False)
            mini_ = Corpus(subtokens)
            n_types = mini_.N
            row_tuple = [m_tokens, n_types]

            # calculate legomena k-vector expression of tokens sample
            if dim is not None:
                lego_k = mini_.k
                lego_k[0] = self.N - mini_.N
                rank1 = len(lego_k)  # frequency of most frequent token
                lego_k = lego_k[:dim]  # only care about n < upto
                pad_width = dim - len(lego_k)
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
        if dim is not None:
            colnames += ["lego_" + str(x) for x in range(dim)]
            colnames += ["rank1"]
        df = pd.DataFrame(TTR, columns=colnames)

        # return
        return df

    def iseries(self, m_tokens: np.ndarray) -> np.ndarray:
        """
        Runs infinite series model to predict N = E(M)
        NOTE: Eqn (8b) in https://arxiv.org/pdf/1901.00521.pdf
        :param m: List-like independent variable m, the number of tokens
        :returns n: Array of dependent variables n, the number of types
        """

        # sum over legomena coefficients k
        m_tokens = np.array(m_tokens).reshape(-1)
        x = m_tokens / self.M
        exponents = range(len(self.k))
        terms = np.array([np.power(1 - x, n) for n in exponents])
        k_0 = np.dot(self.k, terms)
        E_x = np.round(self.N - k_0)

        # return
        return E_x

    # apply new log formula y = f(x) for proportions y,x in [0,1]
    def log_formula(self, x: np.ndarray, epsilon: float = None, dim: int = None):
        """
        Applies the normalized sublinear growth x -> y formula [0,1] -> [0,1]
        n_types = N_z * log_formula(m_tokens / M_z)
        NOTE: Eqns (17.*) in https://arxiv.org/pdf/1901.00521.pdf
        :param x: List-like independent variable x, the proportion of tokens
        :param epsilon: Returns limit value within +/- epsilon of singularity
        :param dim: Desired dimension of k-vector, highest n for which to return n-legomena
        :returns E_x: Expected proportion of types for given proportion of tokens x
        :returns k_frac: Expected proportions of n-legomena for n = 0 to 5
        """

        x = np.array(x).reshape(-1)

        # user options
        epsilon = epsilon or self.epsilon
        dim = dim or self.options.dimension

        # TODO: implement generalized formula
        MAX_DIM = 6
        dim = MAX_DIM

        # initialize return vectors
        _E_x = np.zeros((len(x)))
        _k_frac = np.zeros((len(x), MAX_DIM))

        # two singularities @ x=0 & x=1
        x_iszero = np.abs(x) < epsilon  # singularity @ x=0
        x_isone = np.abs(x - 1.0) < epsilon  # singularity @ x=1
        x_isokay = ~np.logical_or(x_iszero, x_isone)
        x = x[x_isokay]

        # initialize results
        k_frac = np.zeros((len(x), MAX_DIM))
        logx = np.log(x)
        x2, x3, x4, x5, x6 = [x ** i for i in range(2, MAX_DIM + 1)]
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
    def predict(self, m_tokens: np.ndarray):
        """
        Applies the log formula to predict number of types as a function of tokens.
        :param m_tokens: Number of tokens (scalar or list-like)
        :returns E_x: Expected proportion of types for given proportion of tokens x
        :returns k_frac: Expected proportions of n-legomena for n = 0 to 5
        """

        # allow scalar
        return_scalar = np.isscalar(m_tokens)
        m_tokens = np.array(m_tokens).reshape(-1)

        # retrieve fitting parameters
        M_z, N_z = self.params

        # scale down to normalized log formula
        x = m_tokens / M_z
        E_x, k_frac = self.log_formula(x)
        E_m = np.round(N_z * E_x)
        k = np.round(N_z * k_frac)

        # allow scalar
        if return_scalar:
            assert len(E_m) == len(k) == 1
            E_m = E_m.squeeze()
            k = k.squeeze()

        # return scaled up predictions
        return E_m, k

    # find M_z, N_z
    def fit(self, optimize: bool = False) -> LogParams:
        """
        Uses observed whole-corpus hapax:type ratio to infer M_z, N_z optimum sample size.
        :param optimize: (optional) Uses scipy.optimize.curve_fit() to refine initial fit.
        :returns: Logarithmic model parameters, as tuple
        """

        # infer z from h_obs
        h_obs = self.k[1] / self.N  # observed hapax frac
        func = lambda x: 1.0 / np.log(x) + 1.0 / (1.0 - x) - h_obs
        z0 = 0.5
        z = fsolve(func, z0)[0]
        M_z = self.M / z
        N_z = self.N * (z - 1) / z / np.log(z)
        M_z, N_z = int(M_z), int(N_z)
        self.params = (M_z, N_z)

        # minimize MSE on random perturbations of M_z, N_z
        if optimize:
            TTR = self.TTR
            func = lambda m, M_z, N_z: N_z * np.log(m / M_z) * m / M_z / (m / M_z - 1)
            xdata = TTR.m_tokens
            ydata = TTR.n_types
            params, _ = curve_fit(func, xdata, ydata, p0=(M_z, N_z))
            M_z, N_z = int(params[0]), int(params[1])
            self.params = (M_z, N_z)

        # return optimum sample size M_z, N_z
        return self.params


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
        :param pgid: Project Gutenberg book ID
        :returns: Corpus, wrapper class for the "bag of words"
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

        # build corpus from frequency distribution
        df = pd.read_csv(fobj, header=-1, delimiter="\t")
        asdict = {row[1]: row[2] for row in df.itertuples()}
        corpus = Corpus(asdict)

        # return corpus object
        return corpus

    # get SPGC metadata from csv
    def getMeta():
        """Retrieves metadata.csv and returns as dataframe."""

        # read csv
        fname = f"{DATAPATH}/SPGC-metadata-2018-07-18.csv"
        df = pd.read_csv(fname)

        # strip out "sound" entries
        df = df[df.type == "Text"]
        assert df.shape == (56466, 9), f"ERROR: Corrupted SPGC file {fname}."

        # key records by numeric PG ID
        df.id = df.id.str.replace("PG", "")
        df = df.set_index("id")

        # return
        return df
