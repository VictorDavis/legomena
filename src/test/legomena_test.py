
# bloody dependencies
import inspect
from nltk.corpus import gutenberg
import numpy as np
import os
from sklearn.metrics import mean_squared_error as MSE
import unittest

from legomena import Corpus

GRAPHICS_ON = True

# display test name
def print_test_name():
    print("\n"+inspect.stack()[1][3]+"\n", end='', flush=True)

# root mean squared error as a percent of realization
def RMSE_pct(y, y_hat):
    y = np.array(y).reshape(1)
    y_hat = np.array(y_hat).reshape(1)
    return np.sqrt(MSE(y, y_hat)) / y

# unit tests
class LegomenaTest(unittest.TestCase):

    # initialize Legomena class with small corpus & assert simple stats match
    def test_basic(self):
        print_test_name()

        # retrieve corpus from NLTK
        filename = 'blake-poems.txt'
        words = list(gutenberg.words(filename))

        # initialize class
        corpus = Corpus(words)

        # sanity checks
        assert corpus.M == len(words)
        assert corpus.N == len(set(words))
        assert sum(corpus.fdist.values()) == corpus.M
        assert len(corpus.fdist) == corpus.N
        hapaxes = corpus.hapaxes()
        assert corpus.k[0] == 0
        assert corpus.k[1] == len(hapaxes)
        assert sum(corpus.k) == corpus.N

    # model-fitting tests
    def test_models(self):
        print_test_name()

        # retrieve corpus from NLTK
        filename = 'blake-poems.txt'
        words = list(gutenberg.words(filename))

        # initialize class
        corpus = Corpus(words)

        # build TTR curve
        corpus.buildTTRCurve(seed = 42)
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
            heaps_predictions = [ slope * np.log(m) + intercept for m in df.m_tokens ]
            plt.scatter(np.log(df.m_tokens), np.log(df.n_types))
            plt.plot(np.log(df.m_tokens), heaps_predictions, color = "red")
            plt.title("Heap's Model (K,B) = (%f, %f)" % (K,B))
            plt.xlabel("log(tokens)")
            plt.ylabel("log(types)")
            plt.show()

            # normal graph of Heap's Model
            heaps_predictions = corpus.heaps(df.m_tokens)
            plt.scatter(df.m_tokens, df.n_types)
            plt.plot(df.m_tokens, heaps_predictions, color = "red")
            plt.title("Heap's Model (K,B) = (%f, %f)" % (K,B))
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # Infinite Series Model
            iseries_predictions = corpus.iseries(df.m_tokens)
            plt.scatter(df.m_tokens, df.n_types)
            plt.plot(df.m_tokens, iseries_predictions, color = "red")
            plt.title("Infinite Series Model")
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

    # test legomena models
    def test_legomena(self):
        print_test_name()

        # retrieve corpus from NLTK
        filename = 'blake-poems.txt'
        words = list(gutenberg.words(filename))

        # initialize class
        corpus = Corpus(words)

        # calculate transformation matrix A_x
        D = 5
        A_0 = np.array([ [ int(i == 0) for j in range(D)] for i in range(D) ])
        A_1 = np.array([ [ int(i == j) for j in range(D)] for i in range(D) ])
        assert np.array_equal(corpus.transformationMatrix(0.0, D), A_0)
        assert np.array_equal(corpus.transformationMatrix(1.0, D), A_1)
        print( corpus.transformationMatrix(0.5, D) )

        # transform k to k'
        assert np.array_equal(corpus.transform(1), corpus.k) # sample 100% of the corpus and assert k'=k
        k_0 = corpus.transform(0) # k'(0) = [N, 0, 0, 0, ..., 0]
        assert k_0[0] == corpus.N
        assert sum(k_0) == corpus.N

        # build n-legomena curves
        corpus.buildTTRCurve(seed = 42, legomena_upto = 5)
        df = corpus.TTR

        # generate predictions
        m_choices = df.m_tokens           # counts
        x_choices = m_choices / corpus.M  # proportions
        k_matrix  = np.array([ corpus.transform(x) for x in x_choices ])

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # predicted hapaxes
            predictions = k_matrix[:, 1]
            realization = df.lego_1
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color = "red")
            plt.title("Hapax-Token Relation")
            plt.xlabel("tokens")
            plt.ylabel("hapaxes")
            plt.show()

            # predicted hapaxes proportions
            predictions = k_matrix[:, 1] / (corpus.N - k_matrix[:, 0])
            realization = df.lego_1 / df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color = "red")
            plt.title("Hapax-Token Relation (fraction)")
            plt.xlabel("tokens")
            plt.ylabel("hapax fraction")
            plt.show()

    # functions for finding optimum sample size M_z, N_z
    def test_optimization(self):
        print_test_name()

        # retrieve corpus from NLTK
        filename = 'blake-poems.txt'
        words = list(gutenberg.words(filename))

        # initialize class
        corpus = Corpus(words)

        # use observed hapax:type ratio to infer optimum
        # NOTE: no sampling required to parametrize model
        def function_h(x):
            '''Expected hapax:type ratio for proportion x of an optimum sample.'''
            return 1/np.log(x) + 1/(1-x)

        np.random.seed(42)
        h_space  = np.random.uniform(0.25, 0.75, 100)
        EPS = 1e-8
        for h_obs in h_space:
            h_inverse = corpus.inverse_h(h_obs)
            assert abs(function_h(h_inverse) - h_obs) < EPS

        # infer optimum sample size from observed hapax:type ratio
        corpus.fit()
        M_z, N_z = corpus.M_z, corpus.N_z

        # measure E(x), k(x) and compare
        corpus.buildTTRCurve(seed = 42, legomena_upto = 5)
        df = corpus.TTR

        # generate single prediction
        E_m, k = corpus.predict(corpus.M)
        assert E_m  == corpus.N
        assert k[1] == df.lego_1.max()
        assert RMSE_pct(df.lego_2.max(), k[2]) < 0.02

        # generate vector of predictions
        E_m, k = corpus.predict(df.m_tokens)

        # draw pretty pictures
        if GRAPHICS_ON:
            import matplotlib.pyplot as plt

            # predicted hapaxes
            predictions = E_m
            realization = df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color = "red")
            plt.title("Type-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("types")
            plt.show()

            # predicted hapax fraction
            predictions = k[:, 1] / E_m
            realization = df.lego_1 / df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color = "red")
            plt.title("Dis-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("hapax fraction")
            plt.show()

            # predicted dis legomena fraction
            predictions = k[:, 2] / E_m
            realization = df.lego_2 / df.n_types
            plt.scatter(df.m_tokens, realization)
            plt.plot(df.m_tokens, predictions, color = "red")
            plt.title("Dis-Token Relation (Log Formula)")
            plt.xlabel("tokens")
            plt.ylabel("dis legomena fraction")
            plt.show()

def suite():
    suite = unittest.TestSuite()
    suite.addTest(LegomenaTest('test_basic'))
    suite.addTest(LegomenaTest('test_models'))
    suite.addTest(LegomenaTest('test_legomena'))
    suite.addTest(LegomenaTest('test_optimization'))
    return suite


# run the aforementioned suite
if __name__ == '__main__':
    print("\n----------------------------Legomena Tests---------------------------")
    print("----------At the end of each: dot stands for OK, E for ERROR----------")
    runner = unittest.TextTestRunner()
    runner.run(suite())
