
# bloody dependencies
import collections
import math
from nltk.corpus import gutenberg # not actually doing any real nlp, just using corpora
import numpy as np
import os
import pandas as pd
import random
import re

# Project Gutenberg does not allow programmatic "pinging" of site pages, so for the sake of transparency and reproducibility,
# using nltk's built-in corpora and wordstrings with no "manual" preprocessing
def getWordBagFromNLTK(filename):
    words = list(gutenberg.words(filename))
    meta = [ words.index(bracket) for bracket in ['[', ']'] ]
    meta = ' '.join(words[(meta[0]+1):meta[1]])
    return (words, meta)

# # Converts Project Gutenberg text into "bag of words"
# # Reads from manually downloaded file, includes preprocessing
# # To test for "clean" text: set(filter(re.compile('[^a-z0-9\']').match, words))
# def getWordBag(filename):
#     # books = pd.read_csv('data/books.csv', index_col = 0)
#     # book = books.loc[bookid]
#     with open('books/'+filename) as file:
#         lines = file.readlines()
#     lim = [ lines.index(line) for line in lines if 'PROJECT GUTENBERG EBOOK' in line ]
#     start = lim[0]+1
#     end = lim[1]-1
#     text = ''.join(lines[start:end])
#
#     # meta-data from gutenberg header
#     bookid  = [ re.compile('\[EBook #([0-9]+)\]').findall(line)[0] for line in lines[0:start] if '[EBook #' in line ]
#     title   = [ re.compile('Title: (.*)'        ).findall(line)[0] for line in lines[0:start] if 'Title: '  in line ]
#     author  = [ re.compile('Author: (.*)'       ).findall(line)[0] for line in lines[0:start] if 'Author: ' in line ]
#     # release = [ re.compile('Release Date: (.*)').findall(line)[0] for line in lines[0:start] if 'Release Date: ' in line ][0]
#     bookid = bookid[0].strip() if len(bookid) > 0 else ''
#     title  = title[0].strip() if len(title) > 0 else ''
#     author = author[0].strip() if len(author) > 0 else ''
#     if author != '':
#         title += ' by '+author
#
#     # break out paragraphs
#     text = re.sub('\n\n', '<para>', text)
#     text = re.sub('\n', ' ', text)
#     para = text.split('<para>')
#     para = [para for para in para if re.search('[.?!:",-]$', para)] # require paragraphs to end with puncuation
#     # remove all footnotes & meaningless punctuation: colons, semicolons, em-dashes
#     para = [re.sub('\[[0-9]+\]',' ',para) for para in para]
#     para = [re.sub(',|"|“|”|\*|;|--|—|\(|\)|:|\_|~',' ',para) for para in para]
#     text = ' <para> '.join(para)
#
#     # break out sentences (credit: https://regex101.com/r/nG1gU7/27)
#     line = re.split('(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?|!)\s', text)
#     line = [re.sub('\.$',' <stop>', line) for line in line]
#     line = [re.sub('\?$',' <question>', line) for line in line]
#     line = [re.sub('\!$',' <exclamation>', line) for line in line]
#     text = ' '.join(line)
#
#     # break out words
#     word = [word.lower() for word in text.split(' ') if len(word) > 0]
#     text = ' '.join(word)
#
#     # return results
#     return (word, bookid, title)

# Predicts number of types based on input (m) and free parameters (Mz,Nz)
def modelValues(m, Mz, Nz):
    x = m/Mz
    if (math.fabs(x-1) < .01): # discontinuity @ m=M
        types = Nz
        hapax = Nz/2
        dis = Nz/6
        tris = Nz/12
        tetrakis = Nz/20
        pentakis = Nz/30
    else:
        logx = math.log(x)
        types = Nz*logx*x/(x-1)
        hapax = Nz*(x**2 - logx*x - x)/(x-1)**2
        dis = Nz*(x**3 - 2*logx*x**2 - x)/2/(x-1)**3
        tris = Nz*(2*x**4 - 6*logx*x**3 + 3*x**3 - 6*x**2 + x)/6/(x-1)**4
        tetrakis = Nz*(3*x**5 - 12*logx*x**4 + 10*x**4 - 18*x**3 + 6*x**2 - x)/12/(x-1)**5
        pentakis = Nz*(12*x**6 - 60*logx*x**5 + 65*x**5 - 120*x**4 + 60*x**3 - 20*x**2 + 3*x)/60/(x-1)**6
    return (types, hapax, dis, tris, tetrakis, pentakis)

# At what point do the n-legomena most closely resemble a perfect Zipf Distribution?
# When 1/ln(x)+1/(1-x) = H for H the observed proportion of hapaxes
# Note: Because of the nature of this function, Newton's Method is not ideal.
#       Instead, we use a binary search to find f(x)-H = 0, which
#       works nicely since f(x) decreases monotonically from 0 to inf
def solveForX(H):
    # print('H = {}'.format(H))
    last_x = 0
    x = .99
    dx = .5
    timer = 1
    epsilon = .00000001 # 10e-8 -> min ~27 iterations
    while (dx > epsilon) & (timer < 999):
        if (x == 0) | (x == 1): # f(x) can't be evaluated @0,1
            x += epsilon # jigger x slightly
        fx = 1/math.log(x) + 1/(1-x) # f(x)
        if fx > H:# if f(x) > H then x is too low
            if x + dx == last_x:
                dx = dx/2
            last_x = x
            x += dx
        elif fx < H: # if f(x) < H then x is too high
            if (x - dx == last_x):
                dx = dx/2
            while x - dx <= 0: # do not let x go negative
                dx = dx/2
            last_x = x
            x -= dx
        else: # if f(x) = H then x is just right
            # print('Found x in {} iterations'.format(timer))
            return x
        timer += 1
        # print('(x, f(x)) = ({}, {})'.format(x, fx))
    # print('Found x={} in {} iterations.'.format(round(x,4), timer))
    return x


# Type-Token Curve for varying corpus lengths (empirically measured/observed)
def getTTRCurve(decks):
    df = []
    for tokens in decks['tokens'].unique():
        deck  = decks[decks['tokens'] == tokens]
        deck.set_index('repeats')
        types = sum(deck['ntypes'])
        nlegomena = [ deck.loc[n, 'ntypes'] if n in deck.index else 0 for n in range(0, 6) ]
        df.append((tokens, types, nlegomena[1], nlegomena[2], nlegomena[3], nlegomena[4], nlegomena[5]))
    df = pd.DataFrame(df, columns = ['tokens', 'types', 'hapax', 'dis', 'tris', 'tetrakis', 'pentakis'])
    return df

# Calculate RMSE % on predictions
def calcRMSE(df):
    n = df.shape[0]
    N = df['types'].max()
    heaps  = math.sqrt(sum([ (row['types'] - row['types_pred_heaps'] )**2 for index,row in df.iterrows() ])/n)/N
    iseries = math.sqrt(sum([ (row['types'] - row['types_pred_iseries'])**2 for index,row in df.iterrows() ])/n)/N
    model  = math.sqrt(sum([ (row['types'] - row['types_pred_model'] )**2 for index,row in df.iterrows() ])/n)/N
    return [heaps, iseries, model]

# Heap's Law estimate
def predictHeaps(df):
    df = pd.DataFrame(df) # redeclare
    df['logM'] = [ math.log(M) for M in df['tokens']]
    df['logN'] = [ math.log(N) for N in df['types']]
    bestfit = np.polyfit(df['logM'], df['logN'], 1)
    K = math.exp(bestfit[1])
    beta = bestfit[0]
    predictions = [K*m**beta for m in df['tokens']]
    return predictions

# Infinite Series approach
def predictInfSeries(df, M, N, deck):
    deck = deck.sort_values(by = ['repeats'])
    iseries = []
    for tokens in df['tokens']:
        terms = []
        last_term = 1
        for index,row in deck.iterrows():
            if last_term < .01: # series should converge quickly
                break
            k = row['ntypes']
            s = row['repeats']
            x = 1-tokens/M
            if math.fabs(x) <= 1:
                term = k*x**s
                terms.append(term)
                last_term = term
        if len(terms) > 0:
            iseries.append(N - int(sum(terms)))
        else:
            iseries.append(np.nan)
    return iseries

# returns frequency distribution of words within text
def getFreqDist(words):
    fdist = collections.Counter(words)
    df = pd.DataFrame.from_dict(fdist, orient='index')
    df.columns = ['freq']
    df['word'] = fdist.keys()
    df['rank'] = df['freq'].rank(method='first', ascending = False)
    df = df.sort_values(by = ['rank'])
    # "perfect" Zipf's Law prediction (f1 = N)
    N = len(df.index) # number of types
    df['freq_pred_zipf'] = [ int(N / r) for r in df['rank'] ]
    return df

# Summarizes word frequency distribution into layers or "decks" of "s repeats of k word types"
def getDeck(fdist):
    deck = collections.Counter(fdist['freq']) # frequency distribution of frequencies {6:405} = 405 words appear exactly 6 times
    deck = pd.DataFrame.from_dict(deck, orient = 'index')
    deck.columns = ['ntypes']    # k = number of types
    deck['repeats'] = deck.index # s = number of repeats
    M = sum(fdist['freq'])
    N = len(fdist.index)
    deck['tokens'] = M
    deck['types'] = N
    deck['fraction'] = deck['ntypes'] / deck['types']
    deck['types_pred_zipf'] = [ int(N / n / (n + 1)) for n in deck['repeats'] ]
    deck['fraction_pred_zipf'] = deck['types_pred_zipf'] / N
    return deck

# samples corpus at increments and breaks into decks
# Chosen method: randomly sample words WITHOUT replacement
def getDecks(words, nsamples = 250, method = 3):
    random.seed(2701)
    decks = []
    increment = int(len(words) / nsamples)
    for M in range(increment, len(words), increment):
        if method == 1:
            subset = words[0:M] # first M tokens
        elif method == 2:
            start = random.randint(0, len(words)-M)
            subset = words[start:(start+M)] # random block of M tokens
        elif method == 3:
            subset = random.sample(words, M) # M tokens selected at random
        fdist = getFreqDist(subset)
        deck  = getDeck(fdist)
        decks.append(deck)
    fdist = getFreqDist(words)
    deck  = getDeck(fdist)
    decks.append(deck)
    decks = pd.concat(decks)
    return decks

# Aggregate model predictions alongside empirical observations
def getBookStats(decks, ttr):
    M = decks['tokens'].max()
    N = decks['types'].max()
    deck = pd.DataFrame(decks[decks['tokens'] == M])
    ttr['types_pred_iseries'] = predictInfSeries(ttr, M, N, deck)
    ttr['types_pred_heaps'] = predictHeaps(ttr)
    H = deck.loc[1, 'fraction']
    z = solveForX(H)
    Mz = int(M / z)
    Nz = int(N * (z-1)/math.log(z)/z)
    predictions = [ modelValues(m, Mz, Nz) for m in ttr['tokens'] ]
    for i, col in {0: 'types', 1:'hapax', 2:'dis', 3:'tris', 4:'tetrakis', 5:'pentakis'}.items():
        ttr['{}_pred_model'.format(col)]  = [ elem[i] for elem in predictions ]
        ttr['{}_frac'.format(col)] = ttr['{}'.format(col)] / ttr['types']
        ttr['{}_frac_pred'.format(col)] = ttr['{}_pred_model'.format(col)] / ttr['types_pred_model']
    RMSE = calcRMSE(ttr)
    return (M, N, Mz, Nz, RMSE, ttr)

# create ALL book data & write to disk
def allStatsToDisk():
    books = []
    ttrs = [] # type-token ratio growth
    deckslist = []

    # # get all stats from files on disk, if any
    # for filename in os.listdir('books/'):
    #     # book = books.loc[bookid]
    #     words, bookid, title = getWordBag(filename)
    #     print('Processing ID #{}: {}...'.format(bookid, title))
    #     decks = getDecks(words)
    #     ttr   = getTTRCurve(decks)
    #     M, N, Mz, Nz, RMSE, ttr = getBookStats(decks, ttr)
    #     books.append((bookid, title, 'preprocess', M, N, Mz, Nz, RMSE[0], RMSE[1], RMSE[2]))
    #     # for key, val in {'tokens':M, 'types':N, 'Mz':Mz, 'Nz':Nz, 'RMSE_heaps':RMSE[0], 'RMSE_iseries':RMSE[1], 'RMSE_model':RMSE[2]}.items():
    #     #     books.at[bookid, key] = val
    #     decks['bookid'] = bookid
    #     decks['title']  = title
    #     ttr['bookid'] = bookid
    #     ttr['title']  = title
    #     deckslist.append(decks)
    #     ttrs.append(ttr)

    # get all stats from nltk gutenberg corpus
    bookid = 0
    nbooks = len(gutenberg.fileids())
    for filename in gutenberg.fileids(): # ['bible-kjv.txt', 'blake-poems.txt']: #
        bookid += 1
        words, title = getWordBagFromNLTK(filename)
        print('Processing book #{}/{}: {}...'.format(bookid, nbooks, title))
        try:
            # read in cached decks file, if it exists
            decks = pd.read_csv('data/decks.csv')
            decks = decks[decks['bookid'] == bookid]
        except:
            decks = getDecks(words)
        ttr   = getTTRCurve(decks)
        M, N, Mz, Nz, RMSE, ttr = getBookStats(decks, ttr)
        books.append((bookid, title, M, N, Mz, Nz, RMSE[0], RMSE[1], RMSE[2]))
        decks['bookid'] = bookid
        decks['title']  = title
        ttr['bookid'] = bookid
        ttr['title']  = title
        deckslist.append(decks)
        ttrs.append(ttr)

    # aggregate all data into single dataframes
    books = pd.DataFrame(books, columns = ['bookid', 'title', 'tokens', 'types', 'Mz', 'Nz', 'RMSE_heaps', 'RMSE_iseries', 'RMSE_model'])
    books = books.set_index('bookid')
    ttr = pd.concat(ttrs)
    decks = pd.concat(deckslist)

    # write to disk
    books.to_csv('data/books.csv')
    decks.to_csv('data/decks.csv') # cache "raw" data, intensive to recompute
    ttr.to_csv('data/ttr.csv')

allStatsToDisk()
