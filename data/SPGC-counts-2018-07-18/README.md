# The Standard Project Gutenberg Corpus

This folder contains a 17-book sub-selection of the much larger [Standard Project Gutenberg Corpus](https://arxiv.org/abs/1812.08092) created by Font-Clos, F. & Gerlach, M. These are the minimum number of books to run the unit tests for this package, and also corresponds to NLTK's gutenberg corpus, as documented below.

Download the full 55,000-book 'frozen' version of the corpus (SPGC-2018-07-18) as a Zenodo dataset:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2422560.svg)](https://doi.org/10.5281/zenodo.2422560)

## SPGC vs NLTK

The [Natural Language ToolKit](https://www.nltk.org/) is the gold standard of NLP for python. Unfortunately, it contains a rather measly 18-book corpus from [Project Gutenberg](https://www.gutenberg.org/), a rather pathetic sample size for any robust NLP claims.

The [Standard Project Gutenberg Corpus](https://arxiv.org/abs/1812.08092), created by Font-Clos, F. & Gerlach, M. in 2018, consists of two parts:

- [PGCorpus](https://github.com/pgcorpus/gutenberg): a tool for downloading and tokenizing the entire 55,000-book corpus available from [Project Gutenberg](https://www.gutenberg.org/) at any given time.
- [SPGC-2018-07-18](https://doi.org/10.5281/zenodo.2422560): a 'frozen' version of the corpus as it existed at the time of the tool's creation.

Either of these tools may be used to obtain word frequency distributions for _any_ book on [Project Gutenberg](https://www.gutenberg.org/) for NLP research. Because text tokenization is a non-trivial task involving a myriad of preprocessing microdecisions, any researcher potentially introduces bias into her research by attempting to "roll her own" string tokenizer. Therefore, outsourcing this piece to a neutral third party has three benefits:

- **Reproduceability**: Identical word frequency data may be obtained by author and reviewer from a third party.
- **Objectivity**: Preprocessing microdecisions are made by a third party, so the author can't "fine tune" them (consciously or unconsciously) to a particular hypothesis.
- **Robustness**: (SPGC only) Hypotheses may be tested against a 55,000+ book sample.

Here is a lookup reference for the Project Gutenberg IDs of each book in the `nltk.gutenberg` corpus. Because of differing preprocessing microdecisions, the word counts won't be identical between the two sources. Also, PG#11, _Alice in Wonderland_, is missing from SPGC-2018-07-18 for some reason. Curiouser and curiouser!

```
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
    ("whitman-leaves.txt", 1322)
```
