
# bloody dependencies
import inspect
from nltk.corpus import gutenberg
import os
import unittest

from legomena import Legomena

# display test name
def print_test_name():
    print("\n"+inspect.stack()[1][3]+"\n", end='', flush=True)

# unit tests
class LegomenaTest(unittest.TestCase):

    # initialize Legomena class with small corpus & assert simple stats match
    def test_basic(self):
        print_test_name()

        # retrieve Moby Dick from NLTK
        filename = 'blake-poems.txt'
        corpus = list(gutenberg.words(filename))

        # initialize class
        lego = Legomena(corpus)

        # sanity checks
        assert lego.M == len(corpus)
        assert lego.N == len(set(corpus))
        assert sum(lego.fdist.values()) == lego.M
        assert len(lego.fdist) == lego.N
        hapaxes = [ key for key, val in lego.fdist.items() if val == 1]
        assert lego.k[0] == 0
        assert lego.k[1] == len(hapaxes)
        assert sum(lego.k) == lego.N

def suite():
    suite = unittest.TestSuite()
    suite.addTest(LegomenaTest('test_basic'))
    return suite


# run the aforementioned suite
if __name__ == '__main__':
    print("\n----------------------------Legomena Tests---------------------------")
    print("----------At the end of each: dot stands for OK, E for ERROR----------")
    runner = unittest.TextTestRunner()
    runner.run(suite())
