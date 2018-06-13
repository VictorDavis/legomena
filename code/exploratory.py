
# bloody dependencies
import collections
import math
import numpy as np
import pandas as pd
import random
import re

# Converts Project Gutenberg text into "bag of words"
def getBookText(bookid):
    books = pd.read_csv('books.csv', index_col = 0)
    book = books.loc[bookid]
    file = open('books/'+book.filename)
    book_lines = file.readlines()
    lim = [book_lines.index(line) for line in book_lines if 'PROJECT GUTENBERG EBOOK' in line]
    start = lim[0]+1
    end = lim[1]-1
    book_text = ''.join(book_lines[start:end])

    # break out paragraphs
    book_text = re.sub('\n\n', '<para>', book_text)
    book_text = re.sub('\n', ' ', book_text)
    book_para = book_text.split('<para>')
    book_para = [para for para in book_para if re.search('[.?!:",-]$', para)] # require paragraphs to end with puncuation
    # remove all meaningless punctuation: colons, semicolons, em-dashes
    book_para = [re.sub(',|"|;|--|\(|\)|:',' ',para) for para in book_para]
    book_text = ' <para> '.join(book_para)

    # break out sentences (credit: https://regex101.com/r/nG1gU7/27)
    book_line = re.split('(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?|!)\s', book_text)
    book_line = [re.sub('\.$',' <stop>', line) for line in book_line]
    book_line = [re.sub('\?$',' <question>', line) for line in book_line]
    book_line = [re.sub('\!$',' <exclamation>', line) for line in book_line]
    book_text = ' '.join(book_line)

    # break out words
    book_word = [word.lower() for word in book_text.split(' ') if len(word) > 0]
    book_text = ' '.join(book_word)

    # return results
    return book_text

# Predicts number of types based on input (n) and free parameters (M,N,k)
def modelValues(n, M, N, k=1):
    if k != 1:
        M = M * k
        N = N * (k-1)/math.log(k)
    if (n == M):
        types = N
        hapax = int(N/2)
        tris = int(N/6)
        tetrakis = int(N/12)
    else:
        x = n/M
        types = int(-N*math.log(x)*x/(1-x))
        hapax = int(N*x*(x-math.log(x)-1)/(1-x)**2)
    return (types, hapax)

# Predicts number of types based on input (n) and free parameters (M,N,k)
def hapaxValue(n, M, N, k=1):
    if k != 1:
        M = M * k
        N = N * (k-1)/math.log(k)
    if (n == M):
        result = int(N/2)
    else:
        result = int(-N*math.log(n/M)*n/M/(1-n/M))
    return result

# Fit parameter k to TTR data
# TODO: use RMSE instead
def fitModel(M, N, xvals, yvals):
    k = 1
    dk = .5
    area_correct = sum(yvals)
    last_error = 1
    timer = 1
    while (dk > .0001) & (timer < 9999):
        area_prediction = sum([modelValues(n, M, N, k)[0] for n in xvals])
        error = area_prediction - area_correct
        print('k = {}, error = {}, dk = {}'.format(k, error, dk))
        if (last_error > 0) & (error > 0):
            k = k + dk
        if (last_error < 0) & (error < 0):
            k = k - dk
        if (last_error > 0) & (error < 0):
            dk = dk / 2
            k = k + dk
        if (last_error < 0) & (error > 0):
            dk = dk / 2
            k = k - dk
        last_error = error
        timer += 1
    return k

# Type-Token Curve for varying corpus lengths
def getTTRCurve(decks):
    df = []
    for M in decks['tokens'].unique():
        deck  = decks[decks['tokens'] == M]
        types = sum(deck['ntypes'])
        hapax = deck[deck['repeats'] == 1]['ntypes'].values[0]
        dis   = deck[deck['repeats'] == 2]['ntypes'].values[0]
        tris  = deck[deck['repeats'] == 3]['ntypes'].values[0]
        tetrakis = deck[deck['repeats'] == 4]['ntypes'].values[0]
        df.append((M, types, hapax, dis, tris, tetrakis))
    df = pd.DataFrame(df, columns = ['tokens', 'types', 'hapax', 'dis', 'tris', 'tetrakis'])
    return df

# Taylor Series approach
def predictTaylor(df, M, N, deck):
    deck = deck.sort_values(by = ['repeats'])
    taylor = []
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
            taylor.append(N - int(sum(terms)))
        else:
            taylor.append(np.nan)
    return taylor

# returns frequency distribution of words within text
def getFreqDist(words):
    fdist = collections.Counter(words)
    df = pd.DataFrame.from_dict(fdist, orient='index')
    df.columns = ['freq']
    df['word'] = fdist.keys()
    df['length'] = [1 if ('<' in word) else len(word) for word in df['word']]
    df['rank'] = df['freq'].rank(method='first', ascending = False)
    df = df.sort_values(by = ['rank'])
    # Zipf's Law prediction
    N = len(df.index) # number of types
    df['freq_pred_zipf'] = N / df['rank']
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
    deck['types_pred_zipf'] = [int(N / p / (p + 1)) for p in deck['repeats']]
    deck['fraction_pred_zipf'] = deck['types_pred_zipf'] / N
    return deck

# Takes corpus subsets at increments and breaks into decks
def getDecks(words, increment = 1000, method = 3):
    decks = []
    for M in range(increment, len(words), increment):
        if method == 1:
            subset = words[0:M],    # first M tokens
        elif method == 2:
            start = random.randint(0, len(words)-M)
            subset = words[start:(start+M)], # random block of M tokens
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

# create ALL book data & write to disk
def allStatsToDisk():
    books = pd.read_csv('books.csv', index_col = 0)
    ttrs = [] # type-token ratio growth
    for bookid in [2701]: # books.index: #
        book = books.loc[bookid]
        print('Processing {}...'.format(book['title']))
        words = getBookText(bookid).split(' ')
        M = len(words)
        N = len(set(words))
        decks = getDecks(words)
        deck  = decks[decks['tokens'] == M]
        ttr   = getTTRCurve(decks)
        ttr['types_pred_taylor'] = predictTaylor(ttr, M, N, deck)
        k = fitModel(M, N, ttr['tokens'], ttr['types'])
        predictions = [ modelValues(n, M, N, k) for n in ttr['tokens'] ]
        ttr['types_pred_model'] = [ elem[0] for elem in predictions ]
        ttr['hapax_pred_model'] = [ elem[1] for elem in predictions ]
        books.at[bookid, 'tokens'] = M
        books.at[bookid, 'types'] = N
        books.at[bookid, 'k'] = k
        books.at[bookid, 'M'] = int(M * k)
        books.at[bookid, 'N'] = int(-N * math.log(k)*k/(1-k))
        for df in [decks, ttr]:
            df['bookid'] = bookid
            df['title']  = book['title']
            df['author'] = book['author']
        decks.append(decks)
        ttrs.append(ttr)
    ttr   = pd.concat(ttrs)
    books.to_csv('books.csv')
    decks.to_csv('decks.csv')
    ttr.to_csv('ttr.csv')

allStatsToDisk()
