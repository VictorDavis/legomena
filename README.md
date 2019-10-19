# Legomena

Tool for exploring types, tokens, and n-legomena relationships in text. Based on [Davis 2019](https://arxiv.org/pdf/1901.00521.pdf) [1] research paper.

[![Testing Status](https://github.com/VictorDavis/legomena/workflows/tests/badge.svg)](https://github.com/VictorDavis/legomena/actions)

## Installation

```
pip install legomena
```

## Data Sources

This package may be driven by any data source, but the author has tested two: the [Natural Language ToolKit](https://www.nltk.org/) and/or the [Standard Project Gutenberg Corpus](https://arxiv.org/abs/1812.08092). The former being the gold standard of python NLP applications, but having a rather weak 16-book gutenberg corpus. The latter containing the full 40,000+ book gutenberg corpus, already tokenized and counted. NOTE: The overlap of the two datasets do _not_ agree in their exact type/token counts, their methodology differing, but this package takes type/token counts as raw data and is therefore methodology-agnostic.

```
# moby dick from NLTK
import nltk
nltk.download("gutenberg")
from nltk.corpus import gutenberg
words = gutenberg.words("melville-moby_dick")
corpus = Corpus(words)
assert corpus.M, corpus.N == (260819, 19317)

# moby dick from SPGC
# NOTE: download and unzip https://zenodo.org/record/2422561/files/SPGC-counts-2018-07-18.zip into DATA_FOLDER
fname = "%s/SPGC-counts-2018-07-18/PG2701_counts.txt" % DATA_FOLDER
with open(fname) as f:
    df = pd.read_csv(f, delimiter="\t", header=None, names=["word", "freq"])
    f.close()
wfd = {str(row.word): int(row.freq) for row in df.itertuples()}
corpus = Corpus(wfd)
assert corpus.M, corpus.N == (210258, 16402)
```

## Basic Usage:

Demo notebooks may be found [here](https://github.com/VictorDavis/legomena/tree/master/notebooks). Unit tests may be found [here](https://github.com/VictorDavis/legomena/blob/master/test_legomena.py).

```
# basic properties
corpus.tokens  # list of tokens
corpus.types  # list of types
corpus.fdist  # word frequency distribution dataframe
corpus.WFD  # alias for corpus.fdist
corpus.M  # number of tokens
corpus.N  # number of types
corpus.k  # n-legomena vector
corpus.hapax  # number of hapax legomena, alias for corpus.nlegomena(1)
corpus.dis  # number of dis legomena, alias for corpus.nlegomena(2)
corpus.tris  # number of tris legomena, alias for corpus.nlegomena(3)
corpus.tetrakis  # number of tetrakis legomena, alias for corpus.nlegomena(4)
corpus.pentakis  # number of pentakis legomena, alias for corpus.nlegomena(5)

# advanced properties
corpus.options  # tuple of optional settings
corpus.resolution  # number of samples when calculating TTR curve
corpus.dimension  # n-legomena vector length to pre-compute
corpus.seed  # random number seed for sampling TTR data
corpus.TTR  # type-token ratio dataframe

# basic functions
corpus.nlegomena(n:int)  # list of types occurring exactly n times
corpus.sample(m:int)  # samples m tokens from corpus *without replacement*
corpus.sample(x:float)  # samples proportion x of corpus *without replacement*
```

## Type-Token Models

There are a variety of models in the literature predicting number of types as a function of tokens, the most well-known being [Heap's Law](https://en.wikipedia.org/wiki/Heaps%27_law). Here are a few implemented, overlaid by the `Corpus` class.

```
# three models
model = HeapsModel()  # Heap's Law
model = InfSeriesModel(corpus)  # Infinite Series Model [1]
model = LogModel()  # Logarithmic Model [1]

# model fitting
m_tokens = corpus.TTR.m_tokens
n_types = corpus.TTR.n_types
model.fit(m_tokens, n_types)
predictions = model.fit_predict(m_tokens, n_types)

# model parameters
model.params

# model predictions
predictions = model.predict(m_tokens)
```

## Demo App

Check out the [demo app](http://legomena.herokuapp.com/) to explore the type-token and n-legomena counts of a few Project Gutenberg books.
