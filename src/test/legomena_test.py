# bloody dependencies
import inspect
from nltk.corpus import gutenberg
import numpy as np
import os
import pandas as pd
from sklearn.metrics import mean_squared_error as MSE
import unittest

from legomena import Corpus, SPGC

GRAPHICS_ON = True
PGID = 2701  # moby dick

# display test name
def print_test_name():
    print("\n" + inspect.stack()[1][3] + "\n", end="", flush=True)


# standard error
def std_err(y: float, y_hat: float):
    return abs(y_hat - y) / y


# root mean squared error as a percent of total types (normalized)
# NOTE: y = realizations, y_hat = predictions
def RMSE_pct(y, y_hat):
    y = np.array(y)
    y_hat = np.array(y_hat)

    # NOTE: real NRMSE formula does / ymax-ymin but in this application ymin ~= 0
    return np.sqrt(MSE(y, y_hat)) / y.max()


# unit tests
class LegomenaTest(unittest.TestCase):

    # initialize Legomena class with small corpus & assert simple stats match
    def test_basic(self):
        print_test_name()

        # retrieve corpus from NLTK
        filename = "blake-poems.txt"
        words = list(gutenberg.words(filename))

        # initialize class
        corpus = Corpus(words)

        # sanity checks
        assert corpus.M == len(words)
        assert corpus.N == len(set(words))
        assert corpus.fdist.freq.sum() == corpus.M
        assert corpus.fdist.shape[0] == corpus.N
        hapaxes = corpus.hapaxes()
        assert corpus.k[0] == 0
        assert corpus.k[1] == len(hapaxes)
        assert sum(corpus.k) == corpus.N

    # ETL & sample SPGC data to generate legomena counts
    def test_spgc(self):
        print_test_name()

        # retrieve Moby Dick from SPGC
        corpus = SPGC.get(PGID)
        meta = SPGC.getMeta()

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
            print(f"\nLoading SPGC {book.SPGC_id}")
            corpi[f"{book.SPGC_id}_SPGC"] = SPGC.get(book.SPGC_id)
            print(f"\nLoading NLTK {book.NLTK_id}")
            corpi[f"{book.NLTK_id}_NLTK"] = Corpus(list(gutenberg.words(book.NLTK_id)))

        # fit TTR curve for all & compare RMSE
        results = []
        for corpus_name, corpus in corpi.items():
            corpus.buildTTRCurve()

            # fit log model
            corpus.fit()
            predictions, _ = corpus.predict(corpus.TTR.m_tokens)
            realizations = corpus.TTR.n_types
            rmse = RMSE_pct(realizations, predictions)
            print(f"RMSE for {corpus_name} is {rmse}.")

            # fit heaps model
            corpus.fitHeaps()
            predictions = corpus.heaps(corpus.TTR.m_tokens)
            rmse2 = RMSE_pct(realizations, predictions)

            results.append(
                (
                    corpus_name,
                    corpus.M,
                    corpus.N,
                    corpus.M_z,
                    corpus.N_z,
                    rmse,
                    corpus.heaps_K,
                    corpus.heaps_B,
                    rmse2,
                )
            )

        # aggregate & analyze results
        results = pd.DataFrame(results)
        results.columns = [
            "name",
            "M",
            "N",
            "M_z",
            "N_z",
            "RMSE_pct",
            "K",
            "B",
            "RMSE_pct_heaps",
        ]
        results["source"] = [name[-4:] for name in results.name]

        print(results)
        print(
            "SPGC/Heaps avg error:",
            results[results.source == "SPGC"].RMSE_pct_heaps.mean(),
        )
        print(
            "NLTK/Heaps avg error:",
            results[results.source == "NLTK"].RMSE_pct_heaps.mean(),
        )
        print(
            "SPGC/Log   avg error:", results[results.source == "SPGC"].RMSE_pct.mean()
        )
        print(
            "NLTK/Log   avg error:", results[results.source == "NLTK"].RMSE_pct.mean()
        )

    # model-fitting tests
    def test_models(self):
        print_test_name()

        # initialize class
        corpus = SPGC.get(PGID)

        # build TTR curve
        corpus.buildTTRCurve(seed=42)
        df = corpus.TTR

        # fit Heap's Law model to TTR curve
        corpus.fitHeaps()
        K, B = corpus.heaps_K, corpus.heaps_B
        # assert corpus.heaps(1000) == 789

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # log-log graph of Heap's Model
            slope, intercept = B, np.log(K)
            heaps_predictions = [slope * np.log(m) + intercept for m in df.m_tokens]
            plt.scatter(np.log(df.m_tokens), np.log(df.n_types))
            plt.plot(np.log(df.m_tokens), heaps_predictions, color="red")
            plt.title("Heap's Model (K,B) = (%f, %f)" % (K, B))
            plt.xlabel("log(tokens)")
            plt.ylabel("log(types)")
            plt.show()

            # normal graph of Heap's Model
            heaps_predictions = corpus.heaps(df.m_tokens)
            plt.scatter(df.m_tokens, df.n_types)
            plt.plot(df.m_tokens, heaps_predictions, color="red")
            plt.title("Heap's Model (K,B) = (%f, %f)" % (K, B))
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # Infinite Series Model
            iseries_predictions = corpus.iseries(df.m_tokens)
            plt.scatter(df.m_tokens, df.n_types)
            plt.plot(df.m_tokens, iseries_predictions, color="red")
            plt.title("Infinite Series Model")
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

    # test legomena models
    def test_legomena(self):
        print_test_name()

        # initialize class
        corpus = SPGC.get(PGID)

        # calculate transformation matrix A_x
        D = 5
        A_0 = np.concatenate((np.ones(D).reshape(1, D), np.zeros((D - 1, D))), axis=0)
        A_1 = np.identity(D)
        assert np.array_equal(corpus.transformationMatrix(0.0, D), A_0)
        assert np.array_equal(corpus.transformationMatrix(1.0, D), A_1)
        print("A(0.5) = \n", corpus.transformationMatrix(0.5, D))

        # transform k to k'
        # NOTE: conputationally, transform matrix can only handle 1024x1024 dimensions, thus output len(k') <= 1024
        k_0 = corpus.transform(0)  # k'(0) = [N, 0, 0, 0, ..., 0]
        k_1 = corpus.transform(1)  # k'(1) = [0, k1, k2, k3, ...]
        assert np.array_equal(
            k_1, corpus.k
        )  # sample 100% of the corpus and assert k'=k
        assert all(k_0[1:] == 0)
        assert k_0[0] == corpus.N

        # build n-legomena curves
        corpus.buildTTRCurve(seed=42, legomena_upto=5)
        df = corpus.TTR

        # generate predictions
        m_choices = df.m_tokens  # counts
        x_choices = m_choices / corpus.M  # proportions
        k_matrix = corpus.transform(x_choices)
        k_matrix2 = np.array([corpus.transform(x) for x in x_choices])
        assert np.allclose(k_matrix, k_matrix2)

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # predicted hapaxes
            predictions = k_matrix[:, 1]
            realization = df.lego_1
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color="red")
            plt.title("Hapax-Token Relation")
            plt.xlabel("tokens")
            plt.ylabel("hapaxes")
            plt.show()

            # predicted hapaxes proportions
            predictions = k_matrix[:, 1] / (corpus.N - k_matrix[:, 0])
            realization = df.lego_1 / df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color="red")
            plt.title("Hapax-Token Relation (fraction)")
            plt.xlabel("tokens")
            plt.ylabel("hapax fraction")
            plt.show()

    # functions for finding optimum sample size M_z, N_z
    def test_optimization(self):
        print_test_name()

        # retrieve corpus from NLTK
        filename = "melville-moby_dick.txt"
        words = list(gutenberg.words(filename))

        # TODO: Why SPGC model quality suffers relative to NLTK ?

        # initialize class
        corpus = Corpus(words)  # SPGC.get(PGID)

        # infer optimum sample size from observed hapax:type ratio
        corpus.fit()
        M_z, N_z = corpus.M_z, corpus.N_z

        # measure E(x), k(x) and compare
        corpus.buildTTRCurve(seed=42, legomena_upto=5)
        df = corpus.TTR

        # generate single prediction
        E_m, k = corpus.predict(corpus.M)
        assert std_err(corpus.N, E_m) < 0.0001
        assert std_err(corpus.k[1], k[1]) < 0.001
        assert std_err(corpus.k[2], k[2]) < 0.005
        assert std_err(corpus.k[3], k[3]) < 0.05
        assert std_err(corpus.k[4], k[4]) < 0.1
        assert std_err(corpus.k[5], k[5]) < 0.1

        # generate vector of predictions
        E_m, k = corpus.predict(df.m_tokens)

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # predicted hapaxes
            predictions = E_m
            realization = df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color="red")
            plt.title("Type-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # predicted hapax fraction
            predictions = k[:, 1] / E_m
            realization = df.lego_1 / df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color="red")
            plt.title("Dis-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("hapax fraction")
            plt.show()

            # predicted dis legomena fraction
            predictions = k[:, 2] / E_m
            realization = df.lego_2 / df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color="red")
            plt.title("Dis-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("dis legomena fraction")
            plt.show()


def suite():
    suite = unittest.TestSuite()
    suite.addTest(LegomenaTest("test_basic"))
    suite.addTest(LegomenaTest("test_spgc"))
    suite.addTest(LegomenaTest("test_models"))
    suite.addTest(LegomenaTest("test_legomena"))
    suite.addTest(LegomenaTest("test_optimization"))
    return suite


# run the aforementioned suite
if __name__ == "__main__":
    print("\n----------------------------Legomena Tests---------------------------")
    print("----------At the end of each: dot stands for OK, E for ERROR----------")
    runner = unittest.TextTestRunner()
    runner.run(suite())
