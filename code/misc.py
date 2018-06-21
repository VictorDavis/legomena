
# bloody dependencies
import collections
import math
from nltk.corpus import gutenberg
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

# Figure 6: tetrakis legomena
def plotTetrakis():
    mini = pd.read_csv('data/mini.csv')
    data = [
        go.Scatter(
            x = mini['tokens'],
            y = mini['types'],
            name = "Actual",
            mode = "markers"
        ),
        go.Scatter(
            x = mini['tokens'],
            y = mini['predicted'],
            name = "Predicted",
            mode = "lines"
        )
    ]
    layout = go.Layout(
        title = "Moby Dick Tetrakis Mini-Corpus",
        xaxis = dict(title = "M (tokens)"),
        yaxis = dict(title = "N (types)")
    )
    fig = go.Figure(data = data, layout = layout)
    pyo.plot(fig)

def plotWFD(filename = 'melville-moby_dick.txt'):
    words = gutenberg.words(filename)
    fdist = collections.Counter(words)
    df = pd.DataFrame.from_dict(fdist, orient='index')
    df.columns = ['freq']
    df['word'] = fdist.keys()
    df['rank'] = df['freq'].rank(method='first', ascending = False)

    data = [go.Scatter(
        x = df['rank'],
        y = df['freq'],
        mode = 'markers'
    )]
    layout = go.Layout(
        title = "Word Frequency Distributions have Fat Tails",
        xaxis = dict(title = "Rank", type = 'log'),
        yaxis = dict(title = "Frequency", type = 'log')
    )
    fig = go.Figure(data = data, layout = layout)
    pyo.plot(fig)
