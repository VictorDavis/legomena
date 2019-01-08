
# bloody dependencies
import inspect
from nltk.corpus import gutenberg
import numpy as np
import os
import unittest

from legomena import Corpus

GRAPHICS_ON = True

# display test name
def print_test_name():
    print("\n"+inspect.stack()[1][3]+"\n", end='', flush=True)

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
            heaps_predictions = [ slope * x + intercept for x in df.log_m ]
            plt.scatter(df.log_m, df.log_n)
            plt.plot(df.log_m, heaps_predictions, color = "red")
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

def suite():
    suite = unittest.TestSuite()
    suite.addTest(LegomenaTest('test_basic'))
    suite.addTest(LegomenaTest('test_models'))
    return suite


# run the aforementioned suite
if __name__ == '__main__':
    print("\n----------------------------Legomena Tests---------------------------")
    print("----------At the end of each: dot stands for OK, E for ERROR----------")
    runner = unittest.TextTestRunner()
    runner.run(suite())
