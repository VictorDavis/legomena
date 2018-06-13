
# bloody dependencies
import collections
import math
import numpy as np
import pandas as pd
import random
import re

# Converts Project Gutenberg text into "bag of words"
# To test for "clean" text: set(filter(re.compile('[^a-z0-9\']').match, words))
def getWordBag(bookid):
    books = pd.read_csv('data/books.csv', index_col = 0)
    book = books.loc[bookid]
    file = open('books/'+book.filename)
    lines = file.readlines()
    lim = [lines.index(line) for line in lines if 'PROJECT GUTENBERG EBOOK' in line]
    start = lim[0]+1
    end = lim[1]-1
    text = ''.join(lines[start:end])

    # break out paragraphs
    text = re.sub('\n\n', '<para>', text)
    text = re.sub('\n', ' ', text)
    para = text.split('<para>')
    para = [para for para in para if re.search('[.?!:",-]$', para)] # require paragraphs to end with puncuation
    # remove all footnotes & meaningless punctuation: colons, semicolons, em-dashes
    para = [re.sub('\[[0-9]+\]',' ',para) for para in para]
    para = [re.sub(',|"|“|”|\*|;|--|—|\(|\)|:|\_|~',' ',para) for para in para]
    text = ' <para> '.join(para)

    # break out sentences (credit: https://regex101.com/r/nG1gU7/27)
    line = re.split('(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?|!)\s', text)
    line = [re.sub('\.$',' <stop>', line) for line in line]
    line = [re.sub('\?$',' <question>', line) for line in line]
    line = [re.sub('\!$',' <exclamation>', line) for line in line]
    text = ' '.join(line)

    # break out words
    word = [word.lower() for word in text.split(' ') if len(word) > 0]
    text = ' '.join(word)

    # return results
    return word

# Fits hapax_frac data to 1/ln(x)+1/(1-x)
def findOptimum(ttr):
    M = ttr['tokens'].max()
    N = ttr['types'].max()
    Mz = M+1
    dM = M/2
    timer = 1
    realY = ttr['hapax_frac']
    lastError = 1
    while (math.fabs(dM) > 100) & (timer < 999):
        X = ttr['tokens'] / Mz
        predY = [ 1/math.log(x)+1/(1-x) for x in X ]
        totalError = sum((predY - realY)**2)
        if totalError <= lastError:
            while Mz - dM < 0:
                dM = int(dM/2)
        else:
            dM = -int(dM/2)
        Mz = Mz + dM
        lastError = totalError
        print('Mz = {}, dM = {}, Err = {}'.format(Mz, dM, totalError))
        timer += 1
    z = Mz/M
    if (Mz > M): # forecast Nz
        Nz = int(N * math.log(z) * z / (z-1))
    else:
        # lookup nearest approximate value for Nz
        ttr['dM'] = [ math.fabs(m - Mz) for m in ttr['tokens'] ]
        id = ttr['dM'].idxmin()
        Nz = ttr.loc[id]['types']
    return (z, Mz, Nz)

# Predicts number of types based on input (n) and free parameters (Mz,Nz)
def modelValues(n, Mz, Nz):
    x = n/Mz
    if (math.fabs(x-1) < .01): # discontinuity @ n=M
        types = Nz
        hapax = int(Nz/2)
        dis = int(Nz/6)
        tris = int(Nz/12)
        tetrakis = int(Nz/20)
        pentakis = int(Nz/30)
    else:
        logx = math.log(x)
        types = int(Nz*logx*x/(x-1))
        hapax = int(Nz*(x**2 - logx*x - x)/(x-1)**2)
        dis = int(Nz*(x**3 - 2*logx*x**2 - x)/2/(x-1)**3)
        tris = int(Nz*(2*x**4 - 6*logx*x**3 + 3*x**3 - 6*x**2 + x)/6/(1-x)**4)
        tetrakis = int(Nz*(3*x**5 - 12*logx*x**4 + 10*x**4 - 18*x**3 + 6*x**2 - x)/12/(x-1)**5)
        pentakis = int(Nz*(12*x**6 - 60*logx*x**5 + 65*x**5 - 120*x**4 + 60*x**3 - 20*x**2 + 3*x)/60/(x-1)**6)
    return (types, hapax, dis, tris, tetrakis, pentakis)

# # Fit parameter z to TTR data
# # TODO: use RMSE instead
# def fitModel(M, N, xvals, yvals):
#     z = 1
#     dk = .5
#     area_correct = sum(yvals)
#     last_error = 1
#     timer = 1
#     while (dk > .0001) & (timer < 999):
#         area_prediction = sum([modelValues(n, M, N, z)[0] for n in xvals])
#         error = area_prediction - area_correct
#         if error == 0:
#             return z
#         print('z = {}, error = {}, dk = {}'.format(z, error, dk))
#         if (last_error > 0) & (error > 0):
#             z = z + dk
#         if (last_error < 0) & (error < 0):
#             z = z - dk
#         if (last_error > 0) & (error < 0):
#             dk = dk / 2
#             z = z + dk
#         if (last_error < 0) & (error > 0):
#             dk = dk / 2
#             z = z - dk
#         last_error = error
#         timer += 1
#     return z

# At what point do the n-legomena most closely resemble a perfect Zipf Distribution?
# When 1/ln(x)+1/(x-1) = H for H the observed proportion of hapaxes
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


# Type-Token Curve for varying corpus lengths (measured/observed)
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
    taylor = math.sqrt(sum([ (row['types'] - row['types_pred_taylor'])**2 for index,row in df.iterrows() ])/n)/N
    model  = math.sqrt(sum([ (row['types'] - row['types_pred_model'] )**2 for index,row in df.iterrows() ])/n)/N
    return [heaps, taylor, model]

# Calculate distance from the predicted optimum "closure of speech" M ~ N(γ+ln(N))
def calcOffClosure(Mz, Nz):
    minOC = -1
    γ = 0.5772156649
    for N in range(int(.5*Nz), int(1.5*Nz), int(.001*Nz)):
        M = N*(γ+math.log(N))
        OC = (Mz-M)**2+(Nz-N)**2
        if (minOC < 0) | (OC < minOC):
            minOC = OC
    return int(math.sqrt(minOC))

# Heap's Law estimate
def predictHeaps(df):
    df['logM'] = [ math.log(M) for M in df['tokens']]
    df['logN'] = [ math.log(N) for N in df['types']]
    bestfit = np.polyfit(df['logM'], df['logN'], 1)
    K = math.exp(bestfit[1])
    beta = bestfit[0]
    predictions = [K*n**beta for n in df['tokens']]
    return predictions

# # Stepping-stone verification that the card-deck equation can be applied to a mini-corpus
# def minicorpus(words, s = 4, increment = 100):
#     fdist = getFreqDist(words)
#     s_words = fdist[fdist['freq'] == s]['word'].values
#     k = len(set(s_words))
#     ks = k*s
#     words = [word for word in words if word in s_words]
#     df = []
#     for n in range(0, len(words)+1, increment):
#         subset = random.sample(words, n)
#         tokens = len(subset)
#         types  = len(set(subset))
#         predicted = int(k-k*(1-n/ks)**s)
#         df.append((tokens, types, predicted))
#     df.append((ks, k, k))
#     df = pd.DataFrame(df, columns = ['tokens', 'types', 'predicted'])
#     return df

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

# Takes corpus subsets at increments and breaks into decks
# Chosen method: randomly sample words WITHOUT replacement
def getDecks(words, nsamples = 250, method = 3):
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

# create ALL book data & write to disk
def allStatsToDisk():
    books = pd.read_csv('data/books.csv', index_col = 0)
    ttrs = [] # type-token ratio growth
    deckslist = []
    for bookid in books.index: #[16328]: #
        book = books.loc[bookid]
        print('Processing ID #{}: {}...'.format(bookid, book['title']))
        try:
            decks = pd.read_csv('data/decks.csv', index_col = 0)
            decks = pd.DataFrame(decks[decks['bookid'] == bookid])
        except:
            words = getWordBag(bookid)
            decks = getDecks(words)
        M = decks['tokens'].max()
        N = decks['types'].max()
        deck = pd.DataFrame(decks[decks['tokens'] == M])
        try:
            ttr = pd.read_csv('data/ttr.csv', index_col = 0)
            ttr = pd.DataFrame(ttr[ttr['bookid'] == bookid])
        except:
            ttr   = getTTRCurve(decks)
        ttr['types_pred_taylor'] = predictTaylor(ttr, M, N, deck)
        ttr['types_pred_heaps'] = predictHeaps(ttr)
        #z = fitModel(M, N, ttr['tokens'], ttr['types'])
        H = deck.loc[1, 'fraction']
        z = 1/solveForX(H)
        Mz = int(M * z)
        if Mz > M:
            Nz = int(-N * math.log(z)*z/(1-z)) # forecast Nz
        else:
            # lookup nearest approximate value for Nz
            ttr['dM'] = [ math.fabs(m - Mz) for m in ttr['tokens'] ]
            id = ttr['dM'].idxmin()
            Nz = ttr.loc[id]['types']
            ttr = ttr.drop(columns = ['dM'])
        N_pred = modelValues(M, Mz, Nz)[0]
        Nz *= N / N_pred
        predictions = [ modelValues(n, Mz, Nz) for n in ttr['tokens'] ]
        for i, col in {0: 'types', 1:'hapax', 2:'dis', 3:'tris', 4:'tetrakis', 5:'pentakis'}.items():
            ttr['{}_pred_model'.format(col)]  = [ elem[i] for elem in predictions ]
            ttr['{}_frac'.format(col)] = ttr['{}'.format(col)] / ttr['types']
            ttr['{}_frac_pred'.format(col)] = ttr['{}_pred_model'.format(col)] / ttr['types_pred_model']
        RMSE = calcRMSE(ttr)
        # OC = calcOffClosure(Mz, Nz)
        for key, val in {'tokens':M, 'types':N, 'Mz':Mz, 'Nz':Nz, 'RMSE_heaps':RMSE[0], 'RMSE_taylor':RMSE[1], 'RMSE_model':RMSE[2]}.items():
            books.at[bookid, key] = val
        for df in [decks, ttr]:
            df['bookid'] = bookid
            df['title']  = book['title']
            df['author'] = book['author']
        deckslist.append(decks)
        ttrs.append(ttr)
    ttr = pd.concat(ttrs)
    decks = pd.concat(deckslist)
    books.to_csv('data/books.csv')
    decks.to_csv('data/decks.csv')
    ttr.to_csv('data/ttr.csv')

# allStatsToDisk()
