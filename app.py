# bloody dependencies
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
from nltk.corpus import gutenberg
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# internal dependencies
from legomena import Corpus, HeapsModel, InfSeriesModel, LogModel

# model comparison over all books, created by test_spgc_nltk()
books = pd.read_csv("data/books.csv", index_col=0)

# use random seed to ensure cross-figure consistency
SEED = 42


def spgc_read(fileid):

    # open file
    fname = "data/SPGC-counts-2018-07-18/%s" % fileid
    with open(fname) as f:
        df = pd.read_csv(f, delimiter="\t", header=None, names=["word", "freq"])
        f.close()

    # load as dictionary
    wfd = {str(row.word): int(row.freq) for row in df.itertuples()}
    return wfd


# available sources
def sources():
    return [
        {"label": "NLTK - Natural Language ToolKit", "value": "NLTK"},
        {"label": "SPGC - Standard Project Gutenberg Corpus", "value": "SPGC"},
    ]


# Show plots in dashboard
app = dash.Dash(__name__)
app.title = "Legomena Explorer"
server = app.server

app.layout = html.Div(
    [
        html.Div(
            [
                dcc.Dropdown(id="source-selector", options=sources(), value="NLTK"),
                dcc.Dropdown(id="book-selector"),
                dcc.Graph(id="figTTR"),
                dcc.Graph(id="figLego"),
            ],
            style=dict(width="50%", display="inline-block"),
        ),
        html.Div(
            [dcc.Graph(id="figData"), dcc.Graph(id="figRMSE")],
            style=dict(width="50%", display="inline-block"),
        ),
    ]
)


@app.callback(Output("book-selector", "options"), [Input("source-selector", "value")])
def available_books(source):

    # {id:title} by source
    df = books.query("source == @source")
    options = [
        {"label": "%s - %s" % (id, row.title), "value": id} for id, row in df.iterrows()
    ]

    # return
    return options


@app.callback(Output("book-selector", "value"), [Input("book-selector", "options")])
def set_default_book(options):
    return options[0]["value"]


# Fig TTR: Actual vs Predicted types = f(tokens)
@app.callback(
    Output("figTTR", "figure"),
    [Input("source-selector", "value"), Input("book-selector", "value")],
)
def plotTTR(source, fileid):

    # retrieve corpus
    if source == "NLTK":
        words = gutenberg.words(fileid)
        corpus = Corpus(words)
    else:
        wfd = spgc_read(fileid)
        corpus = Corpus(wfd)

    # type-token ratio data
    corpus.seed = SEED
    TTR = corpus.TTR
    m_tokens = TTR.m_tokens
    n_types = TTR.n_types

    # fit Heap's Law model to TTR curve
    hmodel = HeapsModel().fit(m_tokens, n_types)
    TTR["n_types_pred_heaps"] = hmodel.predict(m_tokens)

    # infinite series
    imodel = InfSeriesModel(corpus)
    TTR["n_types_pred_iseries"] = imodel.predict(m_tokens)

    # fit logarithmic model to TTR curve
    lmodel = LogModel().fit(m_tokens, n_types)
    TTR["n_types_pred_log"] = lmodel.predict(m_tokens)
    k_pred = lmodel.predict_k(TTR.m_tokens, dim=corpus.dimension)
    dim = k_pred.shape[1]
    for n in range(dim):
        TTR["lego_%s_pred" % n] = k_pred[:, n]

    # build figure
    data = [
        go.Scatter(
            x=TTR.m_tokens,
            y=TTR[colname],
            name=displayname,
            mode="lines" if "_pred" in colname else "markers",
        )
        for colname, displayname in {
            "n_types": "Types Observed",
            "n_types_pred_heaps": "Types Predicted (Heap's Law)",
            "n_types_pred_iseries": "Types Predicted (Infinite Series)",
            "n_types_pred_log": "Types Predicted (Logarithmic Model)",
            "lego_1": "Hapax Legomena Observed",
            "lego_1_pred": "Hapax Legomena Predicted",
            "lego_2": "Dis Legomena Observed",
            "lego_2_pred": "Dis Predicted",
            "lego_3": "Tris Legomena Observed",
            "lego_3_pred": "Tris Predicted",
            "lego_4": "Tetrakis Legomena Observed",
            "lego_4_pred": "Tetrakis Predicted",
            "lego_5": "Pentakis Legomena Observed",
            "lego_5_pred": "Pentakis Predicted",
        }.items()
    ]
    data.append(
        go.Scatter(
            x=[lmodel.M_z],
            y=[lmodel.N_z],
            name="Optimum Sample",
            marker=dict(
                symbol="circle",
                color="rgba(0,0,0, 0.0)",
                size=42,
                line=dict(color="red", width=2),
            ),
        )
    )
    layout = go.Layout(
        title="Type-Token Ratio & N-Legomena Counts",
        xaxis=dict(title="Corpus Size (Tokens)"),
        yaxis=dict(title="Vocabulary Size (Types)"),
    )
    fig = {"data": data, "layout": layout}
    return fig


# Fig Lego: Hapax & n-legomena proportions as sample size grows
@app.callback(
    Output("figLego", "figure"),
    [Input("source-selector", "value"), Input("book-selector", "value")],
)
def plotLego(source, fileid):

    # retrieve corpus
    if source == "NLTK":
        words = gutenberg.words(fileid)
        corpus = Corpus(words)
    else:
        wfd = spgc_read(fileid)
        corpus = Corpus(wfd)

    # type-token ratio data
    corpus.seed = SEED
    TTR = corpus.TTR
    m_tokens = TTR.m_tokens
    n_types = TTR.n_types

    # fit logarithmic model to TTR curve
    lmodel = LogModel().fit(m_tokens, n_types)
    k_pred = lmodel.predict_k(m_tokens, dim=corpus.dimension, normalize=True)
    dim = k_pred.shape[1]
    for n in range(dim):
        actual = "lego_%s" % n
        predicted = "lego_%s_pred" % n
        TTR[actual] = TTR[actual] / TTR.n_types
        TTR[predicted] = k_pred[:, n]

    # build figure
    data = [
        go.Scatter(
            x=TTR.m_tokens,
            y=TTR[colname],
            name=displayname,
            mode="lines" if "_pred" in colname else "markers",
        )
        for colname, displayname in {
            "lego_1": "Hapax Legomena Observed",
            "lego_1_pred": "Hapax Legomena Predicted",
            "lego_2": "Dis Legomena Observed",
            "lego_2_pred": "Dis Predicted",
            "lego_3": "Tris Legomena Observed",
            "lego_3_pred": "Tris Predicted",
            "lego_4": "Tetrakis Legomena Observed",
            "lego_4_pred": "Tetrakis Predicted",
            "lego_5": "Pentakis Legomena Observed",
            "lego_5_pred": "Pentakis Predicted",
        }.items()
    ]
    data.append(
        go.Scatter(
            x=[lmodel.M_z],
            y=[0.5],
            name="Optimum Sample",
            marker=dict(
                symbol="circle",
                color="rgba(0,0,0, 0.0)",
                size=42,
                line=dict(color="red", width=2),
            ),
        )
    )
    layout = go.Layout(
        title="Type & Legomena Proportions",
        xaxis=dict(title="Corpus Size (Tokens)"),
        yaxis=dict(title="Proportion of Types", tickformat=",.0%"),
    )
    fig = {"data": data, "layout": layout}
    return fig


@app.callback(Output("figData", "figure"), [Input("source-selector", "value")])
def plotData(source):

    # filter to source
    df = books.query("source == @source")
    df = df.sort_values("N_z")

    # M_z = N_z ( H[N_z+1] -1 ) ~ = N_z ( -1 + gamma + log(N_z) )
    df["M_z_pred"] = df.N_z * (-1 + np.euler_gamma + np.log(df.N_z + 1))

    # build figure
    data = [
        go.Scatter(
            x=df.m_tokens, y=df.n_types, name="Actual", text=df.title, mode="markers"
        ),
        go.Scatter(x=df.M_z, y=df.N_z, name="Optimum", text=df.title, mode="markers"),
        go.Scatter(x=df.M_z_pred, y=df.N_z, mode="lines", name="Harmony"),
    ]
    layout = go.Layout(
        title="Actual vs Optimum Type/Token Counts",
        xaxis=dict(title="Corpus Size (Tokens)"),
        yaxis=dict(title="Vocabulary Size (Types)"),
        hovermode="closest",
    )
    fig = {"data": data, "layout": layout}
    return fig


@app.callback(Output("figRMSE", "figure"), [Input("source-selector", "value")])
def plotRMSE(source):
    df = books.query("source == @source")
    data = [
        go.Scatter(x=df.m_tokens, y=df[col], name=label, text=df.title, mode="markers")
        for col, label in {
            "RMSE_heaps": "Heap's Law",
            "RMSE_iseries": "Infinite Series",
            "RMSE_log": "Logarithmic Model",
        }.items()
    ]
    layout = go.Layout(
        title="Root Mean Square Error (% of Types)",
        xaxis=dict(title="Corpus Size (Tokens)"),
        yaxis=dict(title="RMSE %"),
    )
    fig = {"data": data, "layout": layout}
    return fig


# NOTE: to kill "OSError: [Errno 98] Address already in use"
# sudo lsof -t -i tcp:8050 | xargs kill -9

if __name__ == "__main__":
    app.title = "TESTING: " + app.title
    app.run_server(debug=True)
