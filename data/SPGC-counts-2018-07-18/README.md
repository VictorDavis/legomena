# The Standard Project Gutenberg Corpus

This folder contains a 17-book sub-selection of the much larger [Standard Project Gutenberg Corpus](https://arxiv.org/abs/1812.08092) created by Font-Clos, F. & Gerlach, M. These are the minimum number of books to run the unit tests for this package, and also corresponds to NLTK's gutenberg corpus, as documented below.

Download the full 55,000-book 2018-07-18 wordcount snapshot from [Zenodo](https://zenodo.org/record/2422561/files/SPGC-counts-2018-07-18.zip).

It is _not_ necessary to unzip this file in order to use this package. The `SPGC.get(<pgid>)` function will work the same whether `data/SPGC-counts-2018-07-18.zip` exists or `data/SPGC-counts-2018-07-18/PG<pgid>_counts.txt` exists. Therefore, it is _recommended_ that you download this 1.5GB wordcount zipfile before using, adapting, or customizing this code, but not strictly necessary.

## SPGC vs NLTK

The Standard Project Gutenberg Corpus (SPGC) contains a snapshot of _all_ books in the Project Gutenberg library. To download the latest snapshot, check out this repo: https://github.com/pgcorpus/gutenberg. The Natural Language ToolKit (NLTK) sadly only contains 18 books. Here is a lookup reference for the Project Gutenberg IDs of each book in the `nltk.gutenberg` corpus. Because of different preprocessing decisions, the word counts won't be identical between the two sources. For some reason, PG#11, _Alice in Wonderland_, is missing from the SPGC zipfile. Curiouser and curiouser!

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
