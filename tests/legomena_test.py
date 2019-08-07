# bloody dependencies
from nltk.corpus import gutenberg
import numpy as np
import pandas as pd
import unittest

from ..legomena import Corpus, SPGC, HeapsModel, KTransformer, LogModel

GRAPHICS_ON = False
PGID = 2701  # moby dick

# standard error
def std_err(y: float, y_hat: float):
    return abs(y_hat - y) / y


# root mean squared error as a percent of total types (normalized)
# NOTE: y = realizations, y_hat = predictions
def RMSE_pct(y, y_hat):
    y = np.array(y)
    y_hat = np.array(y_hat)

    # NOTE: real NRMSE formula does / ymax-ymin but in this application ymin ~= 0
    RMSE = np.sqrt(np.mean((y - y_hat) ** 2))
    return RMSE / y.max()


# unit tests
class LegomenaTest(unittest.TestCase):

    # initialize Legomena class with small corpus & assert simple stats match
    def test_basic(self):

        # retrieve corpus from NLTK
        filename = "blake-poems.txt"
        words = list(gutenberg.words(filename))

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

    def test_spgc_basic(self):

        # get()
        corpus = SPGC.get(2701)
        corpus2 = SPGC.get("PG2701")
        assert corpus2 == corpus
        assert corpus.M == 210258
        assert corpus.N == 16402
        with self.assertRaises(Exception) as context:
            corpus = SPGC.get(999999)

        # metadata()
        df = SPGC.metadata()
        assert df.shape == (56466, 8)
        df = SPGC.metadata(language="en")
        assert df.shape == (45790, 8)
        df = SPGC.metadata(language="xx")
        assert df.shape == (0, 8)

    # ETL & sample SPGC data to generate legomena counts
    def test_spgc(self):

        # retrieve Moby Dick from SPGC
        corpus = SPGC.get(PGID)
        meta = SPGC.metadata()

        # NLTK-SPGC lookup
        books = pd.DataFrame(
            [
                ("austen-emma.txt", 158),
                ("austen-persuasion.txt", 105),
                ("austen-sense.txt", 161),
                ("bible-kjv.txt", 10),
                ("blake-poems.txt", 574),
                ("bryant-stories.txt", 473),
                ("burgess-busterbrown.txt", 22816),
                # ("carroll-alice.txt", 11),
                ("chesterton-ball.txt", 5265),
                ("chesterton-brown.txt", 223),
                ("chesterton-thursday.txt", 1695),
                ("edgeworth-parents.txt", 3655),
                ("melville-moby_dick.txt", 2701),
                ("milton-paradise.txt", 26),
                ("shakespeare-caesar.txt", 1120),
                ("shakespeare-hamlet.txt", 1122),
                ("shakespeare-macbeth.txt", 1129),
                ("whitman-leaves.txt", 1322),
            ]
        )
        books.columns = ["NLTK_id", "SPGC_id"]

        # build corpi from each source
        corpi = {}
        for book in books.itertuples():
            corpi[f"{book.SPGC_id}_SPGC"] = SPGC.get(book.SPGC_id)
            words = list(gutenberg.words(book.NLTK_id))
            corpi[f"{book.NLTK_id}_NLTK"] = Corpus(words)

        # fit TTR curve for all & compare RMSE
        results = []
        for corpus_name, corpus in corpi.items():
            corpus.seed = 42
            m_tokens = corpus.TTR.m_tokens.values
            n_types = corpus.TTR.n_types.values
            realizations = n_types

            # log model
            predictions = LogModel().fit_predict(m_tokens, n_types)
            rmse = RMSE_pct(realizations, predictions)

            # heaps model
            predictions = HeapsModel().fit_predict(m_tokens, n_types)
            rmse2 = RMSE_pct(realizations, predictions)

            results.append((corpus_name, rmse, rmse2))

        # aggregate & analyze results
        results = pd.DataFrame(results)
        results.columns = ["name", "RMSE_pct", "RMSE_pct_heaps"]
        results["source"] = [name[-4:] for name in results.name]

        # assert log model outperforms heaps
        assert results.RMSE_pct.max() < 0.01
        assert results.RMSE_pct_heaps.min() > 0.01

    # model-fitting tests
    def test_models(self):

        # initialize class
        corpus = SPGC.get(PGID)

        # build TTR curve
        corpus.seed = 42
        TTR = corpus.TTR
        m_tokens = TTR.m_tokens.values
        n_types = TTR.n_types.values

        # fit Heap's Law model to TTR curve
        hmodel = HeapsModel().fit(m_tokens, n_types)
        predictions_heaps = hmodel.predict(m_tokens)
        assert hmodel.predict(1000) == 764

        # infinite series
        predictions_iseries = corpus.iseries(TTR.m_tokens)
        assert corpus.iseries(1000) == 513

        # fit logarithmic model to TTR curve
        lmodel = LogModel().fit(m_tokens, n_types)
        predictions_log = lmodel.predict(m_tokens)
        assert lmodel.predict(1000) == 519

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
            plt.title(f"Logarithmic Model (M_z, N_z) = (%s, %s)" % lmodel.params)
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

    # test legomena models
    def test_legomena(self):

        # initialize class
        corpus = SPGC.get(PGID)

        # calculate transformation matrix A_x
        D = 5
        A_0 = np.concatenate((np.ones(D).reshape(1, D), np.zeros((D - 1, D))), axis=0)
        A_1 = np.identity(D)
        assert np.array_equal(KTransformer.A_x(0.0, D), A_0)
        assert np.array_equal(KTransformer.A_x(1.0, D), A_1)
        A_half = KTransformer.A_x(0.5, D)

        # transform k to k'
        # NOTE: computationally, transform matrix can only handle 1024x1024 dimensions, thus output len(k') <= 1024
        k_0 = KTransformer.transform(corpus.k, 0)  # k'(0) = [N, 0, 0, 0, ..., 0]
        k_1 = KTransformer.transform(corpus.k, 1)  # k'(1) = [0, k1, k2, k3, ...]
        assert np.array_equal(
            k_1, corpus.k
        )  # sample 100% of the corpus and assert k'=k
        assert all(k_0[1:] == 0)
        assert k_0[0] == corpus.N

        # build n-legomena curves
        corpus.dimension = 5
        corpus.seed = 42
        TTR = corpus.TTR

        # generate predictions
        m_choices = TTR.m_tokens  # counts
        x_choices = m_choices / corpus.M  # proportions
        k_matrix = KTransformer.transform(corpus.k, x_choices)

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # predicted hapaxes
            predictions = k_matrix[:, 1]
            realization = TTR.lego_1
            plt.scatter(m_tokens, realization)
            plt.plot(m_tokens, predictions, color="red")
            plt.title("Hapax-Token Relation")
            plt.xlabel("tokens")
            plt.ylabel("hapaxes")
            plt.show()

            # predicted hapaxes proportions
            predictions = k_matrix[:, 1] / (corpus.N - k_matrix[:, 0])
            realization = TTR.lego_1 / TTR.n_types
            plt.scatter(m_tokens, realization)
            plt.plot(m_tokens, predictions, color="red")
            plt.title("Hapax-Token Relation (fraction)")
            plt.xlabel("tokens")
            plt.ylabel("hapax fraction")
            plt.show()

    # functions for finding optimum sample size M_z, N_z
    def test_optimization(self):

        # retrieve corpus from NLTK
        filename = "melville-moby_dick.txt"
        words = list(gutenberg.words(filename))

        # TODO: Why SPGC model quality suffers relative to NLTK ?

        # initialize class
        corpus = Corpus(words)  # SPGC.get(PGID)
        corpus.dimension = 5
        corpus.seed = 42
        TTR = corpus.TTR

        # infer optimum sample size from observed hapax:type ratio
        hapax = corpus.k[1]
        model = LogModel().fit_naive(corpus.M, corpus.N, hapax)
        m_tokens, n_types = TTR.m_tokens, TTR.n_types

        # generate single prediction
        E_m = model.predict(corpus.M)
        k = model.predict_k(corpus.M)
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
            k = model.predict_k(m_tokens)
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
