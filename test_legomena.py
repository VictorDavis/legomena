# bloody dependencies
import nltk

nltk.download("gutenberg")
from nltk.corpus import gutenberg
import numpy as np
import pandas as pd
import unittest

# classes to test
from legomena import Corpus, HeapsModel, LogModel, InfSeriesModel

# globals
GRAPHICS_ON = True
SEED = 42

# standard error
def std_err(y, y_hat):
    return np.abs(y_hat - y) / y


# root mean squared error as a percent of total types (normalized)
# NOTE: y = realizations, y_hat = predictions
def RMSE_pct(y, y_hat):
    y = np.array(y)
    y_hat = np.array(y_hat)

    # NOTE: real NRMSE formula does / ymax-ymin but in this application ymin ~= 0
    RMSE = np.sqrt(np.mean((y - y_hat) ** 2))
    return RMSE / y.max()


def spgc_read(fileid):

    # open file
    fname = "data/SPGC-counts-2018-07-18/%s" % fileid
    with open(fname) as f:
        df = pd.read_csv(f, delimiter="\t", header=None, names=["word", "freq"])
        f.close()

    # load as dictionary
    wfd = {str(row.word): int(row.freq) for row in df.itertuples()}
    return wfd


# unit tests
class LegomenaTest(unittest.TestCase):

    # initialize Legomena class with small corpus & assert simple stats match
    def test_basic(self):

        # retrieve corpus from NLTK
        filename = "blake-poems.txt"
        words = gutenberg.words(filename)

        # initialize class
        corpus = Corpus(words)

        # sanity checks
        assert set(corpus.tokens) == set(words)
        assert corpus.M == len(words)
        assert corpus.N == len(set(words))
        assert corpus.fdist.freq.sum() == corpus.M
        assert corpus.fdist.shape[0] == corpus.N
        assert corpus.k[0] == 0
        assert corpus.k[1] == len(corpus.hapax)
        assert corpus.k[2] == len(corpus.dis)
        assert corpus.k[3] == len(corpus.tris)
        assert corpus.k[4] == len(corpus.tetrakis)
        assert corpus.k[5] == len(corpus.pentakis)
        assert corpus.k[42] == len(corpus.nlegomena(42))
        assert sum(corpus.k) == corpus.N
        assert len(corpus.types) == corpus.N
        assert corpus.WFD.equals(corpus.fdist)
        assert corpus.as_datarow(7) == (8354, 1820, 0, 1009, 293, 138, 72, 57, 41)

        # sample()
        assert corpus.sample(99).M == 99
        assert corpus.sample(x=0).M == 0
        assert corpus.sample(x=1).M == corpus.M
        with self.assertRaises(Exception) as context:
            corpus.sample(x=1.5)  # can't over sample with replace=False

    # compare models on SPGC and NLTK data
    def test_spgc_nltk(self):

        # NLTK-SPGC lookup
        books = pd.DataFrame(
            [
                ("austen-emma.txt", "PG158_counts.txt"),
                ("austen-persuasion.txt", "PG105_counts.txt"),
                ("austen-sense.txt", "PG161_counts.txt"),
                ("bible-kjv.txt", "PG10_counts.txt"),
                ("blake-poems.txt", "PG574_counts.txt"),
                ("bryant-stories.txt", "PG473_counts.txt"),
                ("burgess-busterbrown.txt", "PG22816_counts.txt"),
                # ("carroll-alice.txt", "PG11_counts.txt"),  # SPGC missing PG11 :(
                ("chesterton-ball.txt", "PG5265_counts.txt"),
                ("chesterton-brown.txt", "PG223_counts.txt"),
                ("chesterton-thursday.txt", "PG1695_counts.txt"),
                ("edgeworth-parents.txt", "PG3655_counts.txt"),
                ("melville-moby_dick.txt", "PG2701_counts.txt"),
                ("milton-paradise.txt", "PG26_counts.txt"),
                ("shakespeare-caesar.txt", "PG1120_counts.txt"),
                ("shakespeare-hamlet.txt", "PG1122_counts.txt"),
                ("shakespeare-macbeth.txt", "PG1129_counts.txt"),
                ("whitman-leaves.txt", "PG1322_counts.txt"),
            ],
            columns=["nltk_id", "spgc_id"],
        )

        # build corpi from each source
        corpi = {}
        for book in books.itertuples():
            title = book.nltk_id.split(".")[0]
            wfd = spgc_read(book.spgc_id)
            corpus = Corpus(wfd)
            corpi[book.spgc_id] = ("SPGC", title, corpus)
            words = gutenberg.words(book.nltk_id)
            corpus = Corpus(words)
            corpi[book.nltk_id] = ("NLTK", title, corpus)

        # fit TTR curve for all & compare RMSE
        results = []
        for corpus_id, (source, title, corpus) in corpi.items():
            corpus.seed = SEED
            TTR = corpus.TTR
            m_tokens = TTR.m_tokens.values
            n_types = TTR.n_types.values

            # log model
            model = LogModel()
            predictions = model.fit_predict(m_tokens, n_types)
            rmse_log = RMSE_pct(n_types, predictions)

            # iseries model
            predictions = InfSeriesModel(corpus).predict(m_tokens)
            rmse_iseries = RMSE_pct(n_types, predictions)

            # heaps model
            predictions = HeapsModel().fit_predict(m_tokens, n_types)
            rmse_heaps = RMSE_pct(n_types, predictions)

            # append results
            results.append(
                (
                    corpus_id,
                    title,
                    source,
                    corpus.M,
                    corpus.N,
                    model.M_z,
                    model.N_z,
                    rmse_log,
                    rmse_iseries,
                    rmse_heaps,
                )
            )

        # aggregate & analyze results
        results = pd.DataFrame(
            results,
            columns=[
                "id",
                "title",
                "source",
                "m_tokens",
                "n_types",
                "M_z",
                "N_z",
                "RMSE_log",
                "RMSE_iseries",
                "RMSE_heaps",
            ],
        ).set_index("id")

        # save results to data/books.csv
        results.to_csv("data/books.csv")

        # assert log model outperforms heaps
        assert results.RMSE_log.max() < 0.01
        assert results.RMSE_heaps.min() > 0.008

    # model-fitting tests
    def test_models(self):

        # initialize class
        wfd = spgc_read("PG2701_counts.txt")  # moby dick
        corpus = Corpus(wfd)

        # build TTR curve
        corpus.seed = SEED
        TTR = corpus.TTR
        m_tokens = TTR.m_tokens.values
        n_types = TTR.n_types.values

        # fit Heap's Law model to TTR curve
        hmodel = HeapsModel().fit(m_tokens, n_types)
        predictions_heaps = hmodel.predict(m_tokens)
        H = hmodel.predict(1000)

        # infinite series
        imodel = InfSeriesModel(corpus)
        predictions_iseries = imodel.predict(m_tokens)
        I = imodel.predict(1000)

        # fit logarithmic model to TTR curve
        lmodel = LogModel().fit(m_tokens, n_types)
        predictions_log = lmodel.predict(m_tokens)
        L = lmodel.predict(1000)

        # explicit check
        assert (H, I, L) == (756, 513, 515)

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # log-log graph of Heap's Model
            plt.scatter(m_tokens, n_types)
            plt.plot(m_tokens, predictions_heaps, color="red")
            plt.title("Heap's Model (K,B) = (%0.4f, %0.4f)" % hmodel.params)
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("log(tokens)")
            plt.ylabel("log(types)")
            plt.show()

            # normal graph of Heap's Model
            plt.scatter(m_tokens, n_types)
            plt.plot(m_tokens, predictions_heaps, color="red")
            plt.title("Heap's Model (K,B) = (%0.4f, %0.4f)" % hmodel.params)
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # Infinite Series Model
            plt.scatter(m_tokens, n_types)
            plt.plot(m_tokens, predictions_iseries, color="red")
            plt.title("Infinite Series Model")
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # Logarithmic Model
            plt.scatter(m_tokens, n_types)
            plt.plot(m_tokens, predictions_log, color="red")
            plt.title("Logarithmic Model (M_z, N_z) = (%s, %s)" % lmodel.params)
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

    # functions for finding optimum sample size M_z, N_z
    def test_optimization(self):

        # retrieve corpus from NLTK
        filename = "melville-moby_dick.txt"
        words = gutenberg.words(filename)

        # TODO: Why SPGC model quality suffers relative to NLTK ?

        # initialize class
        corpus = Corpus(words)
        dim = 6
        corpus.dimension = dim
        corpus.seed = SEED
        TTR = corpus.TTR

        # infer optimum sample size from observed hapax:type ratio
        hapax = corpus.k[1]
        model = LogModel().fit_naive(corpus.M, corpus.N, hapax)
        m_tokens, n_types = TTR.m_tokens, TTR.n_types

        # generate single prediction
        E_m = model.predict(corpus.M)
        k = model.predict_k(corpus.M, dim)
        assert std_err(corpus.N, E_m) < 0.0001
        assert std_err(corpus.k[1], k[1]) < 0.001
        assert std_err(corpus.k[2], k[2]) < 0.005
        assert std_err(corpus.k[3], k[3]) < 0.05
        assert std_err(corpus.k[4], k[4]) < 0.1
        assert std_err(corpus.k[5], k[5]) < 0.1

        # optimized predictions: worse fit at m=M, better fit overall
        E_m = model.predict(m_tokens)
        RMSE_before = RMSE_pct(n_types, E_m)
        E_m = model.fit_predict(m_tokens, n_types)
        RMSE_after = RMSE_pct(n_types, E_m)
        assert RMSE_after < RMSE_before

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # predicted hapaxes
            predictions = E_m
            realization = n_types
            plt.scatter(m_tokens, realization)
            plt.plot(m_tokens, predictions, color="red")
            plt.title("Type-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # predicted hapax fraction
            k = model.predict_k(m_tokens, dim)
            predictions = k[:, 1] / E_m
            realization = TTR.lego_1 / n_types
            plt.scatter(m_tokens, realization)
            plt.plot(m_tokens, predictions, color="red")
            plt.title("Hapax-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("hapax fraction")
            plt.show()

            # predicted dis legomena fraction
            predictions = k[:, 2] / E_m
            realization = TTR.lego_2 / n_types
            plt.scatter(m_tokens, realization)
            plt.plot(m_tokens, predictions, color="red")
            plt.title("Dis-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("dis legomena fraction")
            plt.show()

    # enforce random seed behavior
    def test_random_seed(self):

        # retrieve corpus from NLTK
        filename = "blake-poems.txt"
        words = gutenberg.words(filename)

        # sampling process works
        corpus = Corpus(words)
        corpus.seed = SEED
        before = corpus.TTR.n_types
        corpus = Corpus(words)
        corpus.seed = SEED
        after = corpus.TTR.n_types
        pd.testing.assert_series_equal(before, after)

        # explicit checks
        assert list(corpus.TTR.n_types[:9].values) == [
            73,
            115,
            176,
            215,
            253,
            285,
            328,
            359,
            381,
        ]
        assert corpus.sample(9).tokens == [
            '"',
            "O",
            "Tongue",
            "a",
            "fear",
            "go",
            "of",
            "the",
            "the",
        ]

        # model fitting works
        TTR = corpus.TTR
        model = LogModel()
        model.fit(TTR.m_tokens, TTR.n_types)
        before = model.params
        model.fit(TTR.m_tokens, TTR.n_types)
        after = model.params
        assert before == after

    # generalized formula matches explicit formula
    def test_generalized(self):

        # apply both formulas to random selection of x values
        model = LogModel()
        np.random.seed(SEED)
        x = 2 * np.random.rand(99)

        # explicit formulas
        k0 = model.formula_0(x)
        k1 = model.formula_1(x)
        k2 = model.formula_2(x)
        k3 = model.formula_3(x)
        k4 = model.formula_4(x)
        k5 = model.formula_5(x)

        # generalized formula
        kn0 = model.formula_n(0, x)
        kn1 = model.formula_n(1, x)
        kn2 = model.formula_n(2, x)
        kn3 = model.formula_n(3, x)
        kn4 = model.formula_n(4, x)
        kn5 = model.formula_n(5, x)

        # all "close"
        assert RMSE_pct(k0, kn0) < 1e-15
        assert RMSE_pct(k1, kn1) < 1e-12
        assert RMSE_pct(k2, kn2) < 1e-10
        assert RMSE_pct(k3, kn3) < 1e-08
        assert RMSE_pct(k4, kn4) < 1e-05
        assert RMSE_pct(k5, kn5) < 1e-03

    # generalized formula contains no singularities
    def test_singularities(self):

        # singularity @ x=0
        model = LogModel()

        # explicit formula singularity @x=0
        x = np.array([0.0])
        assert np.isnan(model.formula_0(x)[0])
        assert np.isnan(model.formula_1(x)[0])
        assert np.isnan(model.formula_2(x)[0])
        assert np.isnan(model.formula_3(x)[0])
        assert np.isnan(model.formula_4(x)[0])
        assert np.isnan(model.formula_5(x)[0])

        # generalized formula returns correct answer
        assert model.formula_n(0, x)[0] == 1.0
        assert model.formula_n(1, x)[0] == 0.0
        assert model.formula_n(2, x)[0] == 0.0
        assert model.formula_n(3, x)[0] == 0.0
        assert model.formula_n(4, x)[0] == 0.0
        assert model.formula_n(5, x)[0] == 0.0

        # singularity @ x=1
        eps = 1e-6
        x = np.array([1 - eps, 1 + eps])

        # limiting values as eps -> 0, x -> 1
        k0_ = 0.0
        k1_ = 1 / 1 - 1 / 2  # 1/2
        k2_ = 1 / 2 - 1 / 3  # 1/6
        k3_ = 1 / 3 - 1 / 4  # 1/12
        k4_ = 1 / 4 - 1 / 5  # 1/20
        k5_ = 1 / 5 - 1 / 6  # 1/30

        # explicit formulas
        explicit_delta0 = abs(model.formula_0(x) - k0_)
        explicit_delta1 = abs(model.formula_1(x) - k1_)
        explicit_delta2 = abs(model.formula_2(x) - k2_)
        explicit_delta3 = abs(model.formula_3(x) - k3_)
        explicit_delta4 = abs(model.formula_4(x) - k4_)
        explicit_delta5 = abs(model.formula_5(x) - k5_)

        # generalized formula
        generalized_delta0 = abs(model.formula_n(0, x) - k0_)
        generalized_delta1 = abs(model.formula_n(1, x) - k1_)
        generalized_delta2 = abs(model.formula_n(2, x) - k2_)
        generalized_delta3 = abs(model.formula_n(3, x) - k3_)
        generalized_delta4 = abs(model.formula_n(4, x) - k4_)
        generalized_delta5 = abs(model.formula_n(5, x) - k5_)

        # explicit formula diverges from limiting value
        assert all(explicit_delta1 > eps)
        assert all(explicit_delta2 > eps * 1e5)
        assert all(explicit_delta3 > eps * 1e13)
        assert all(explicit_delta4 > eps * 1e18)
        assert all(explicit_delta5 > eps * 1e25)

        # generalized formula tends toward limiting value
        assert all(generalized_delta0 < eps)
        assert all(generalized_delta1 < eps)
        assert all(generalized_delta2 < eps)
        assert all(generalized_delta3 < eps)
        assert all(generalized_delta4 < eps)
        assert all(generalized_delta5 < eps)

        # explicit formula singularity @x=1
        x = np.array([1.0])
        assert np.isnan(model.formula_0(x)[0])
        assert np.isnan(model.formula_1(x)[0])
        assert np.isnan(model.formula_2(x)[0])
        assert np.isnan(model.formula_3(x)[0])
        assert np.isnan(model.formula_4(x)[0])
        assert np.isnan(model.formula_5(x)[0])

        # generalized formula returns correct answer
        assert model.formula_n(0, x)[0] == k0_
        assert model.formula_n(1, x)[0] == k1_
        assert model.formula_n(2, x)[0] == k2_
        assert model.formula_n(3, x)[0] == k3_
        assert model.formula_n(4, x)[0] == k4_
        assert model.formula_n(5, x)[0] == k5_

    # check model k_n predictions for n > 5
    def test_high_dimensions(self):

        # retrieve corpus from NLTK
        filename = "melville-moby_dick.txt"
        words = gutenberg.words(filename)

        # initialize class
        corpus = Corpus(words)
        dim = 64
        corpus.dimension = dim
        corpus.seed = SEED
        TTR = corpus.TTR

        # k vector calculated correctly
        raw = corpus.k[:dim]
        agg = TTR.tail(1).values[0, 2:]
        assert all(agg == raw)

        # create & fit model
        m_tokens, n_types = TTR.m_tokens, TTR.n_types
        model = LogModel().fit(m_tokens, n_types)

        # predict E(M) & k_n(M)
        model.dimension = dim
        E_m = model.predict(corpus.M)
        k = model.predict_k(corpus.M, dim)
        assert len(k) == dim
        err = std_err(corpus.k[:dim], k)

        # visualize error
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            plt.plot(err)
            plt.title("Error Analysis")
            plt.xlabel("n = k-vector index")
            plt.ylabel("std error")
            plt.show()
