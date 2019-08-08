# bloody dependencies
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objects as go

# internal dependencies
from legomena import SPGC, HeapsModel, InfSeriesModel, LogModel

# remember to set to False
TESTING = True

# all have to be read into memory?!
books = pd.read_csv("data/books.csv", index_col=0)

# available books
def inventory():
    """Return all books in the Standard Project Gutenberg Corpus"""

    # get metadata
    df = SPGC.metadata()
    if TESTING:
        df = df.head()
    options = [
        {"label": f"{id} - {book.title}", "value": id} for id, book in df.iterrows()
    ]
    return options


# Fig Data: Real vs. Optimum Type/Token scatterplot
def graphData():
    data = [
        go.Scatter(
            x=books["tokens"],
            y=books["types"],
            name="Actual",
            text=books["title"],
            mode="markers",
        ),
        go.Scatter(
            x=books["Mz"],
            y=books["Nz"],
            name="Optimum",
            text=books["title"],
            mode="markers",
        ),
    ]
    layout = go.Layout(
        title="Actual vs Optimum Type/Token Values",
        xaxis=dict(title="Corpus Size (Tokens)"),
        yaxis=dict(title="Vocabulary Size (Types)"),
    )
    fig = {"data": data, "layout": layout}
    return fig


# Fig RMSE: RMSE by corpus size
def graphRMSE():
    data = [
        go.Scatter(
            x=books["tokens"],
            y=books[col],
            name=label,
            text=books["title"],
            mode="markers",
        )
        for col, label in {
            "RMSE_heaps": "Heap's Law",
            "RMSE_iseries": "Infinite Series",
            "RMSE_model": "Logarithmic Model",
        }.items()
    ]
    layout = go.Layout(
        title="Root Mean Square Error (% of Types)",
        xaxis=dict(title="Corpus Size (Tokens)"),
        yaxis=dict(title="RMSE %"),
    )
    fig = {"data": data, "layout": layout}
    return fig


# Show plots in dashboard
app = dash.Dash(__name__)
app.title = "Legomena Explorer"
server = app.server

app.layout = html.Div(
    [
        html.Div(
            [
                dcc.Dropdown(id="book-selector", options=inventory(), value="PG1"),
                dcc.Graph(id="figTTR"),  # TTR Curve
                dcc.Graph(id="figLego"),  # Hapax & n-legomena proportions
            ],
            style=dict(width="50%", display="inline-block"),
        ),
        html.Div(
            [
                dcc.Graph(
                    id="figData", figure=graphData()
                ),  # actual vs. optimum types by corpus size
                dcc.Graph(id="figRMSE", figure=graphRMSE()),  # RMSE by corpus size
            ],
            style=dict(width="50%", display="inline-block"),
        ),
    ]
)

# Fig TTR: Actual vs Predicted types = f(tokens)
@app.callback(Output("figTTR", "figure"), [Input("book-selector", "value")])
def plotTTR(bookid):

    # retrieve corpus
    corpus = SPGC.get(bookid)
    TTR = corpus.TTR
    m_tokens = TTR.m_tokens.values
    n_types = TTR.n_types.values

    # fit Heap's Law model to TTR curve
    hmodel = HeapsModel().fit(m_tokens, n_types)
    TTR["n_types_pred_heaps"] = hmodel.predict(m_tokens)

    # infinite series
    imodel = InfSeriesModel(corpus)
    TTR["n_types_pred_iseries"] = imodel.predict(m_tokens)

    # fit logarithmic model to TTR curve
    lmodel = LogModel().fit(m_tokens, n_types)
    TTR["n_types_pred_log"] = lmodel.predict(m_tokens)
    k_pred = lmodel.predict_k(TTR.m_tokens)
    dim = k_pred.shape[1]
    for n in range(dim):
        TTR[f"lego_{n}_pred"] = k_pred[:, n]

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
    layout = go.Layout(
        title="Type-Token Ratio & N-Legomena Counts",
        xaxis=dict(title="Tokens"),
        yaxis=dict(title="Types"),
        shapes=[
            go.layout.Shape(
                type="line",
                x0=lmodel.M_z,
                y0=0,
                x1=lmodel.M_z,
                y1=1,
                yref="paper",
                line=dict(color="red", width=3),
            ),
            go.layout.Shape(
                type="line",
                x0=0,
                y0=lmodel.N_z,
                x1=1,
                y1=lmodel.N_z,
                xref="paper",
                line=dict(color="red", width=3),
            ),
        ],
    )
    fig = {"data": data, "layout": layout}
    return fig


# Fig Lego: Hapax & n-legomena proportions as sample size grows
@app.callback(Output("figLego", "figure"), [Input("book-selector", "value")])
def plotLego(bookid):

    # retrieve corpus
    corpus = SPGC.get(bookid)
    TTR = corpus.TTR
    m_tokens = TTR.m_tokens.values
    n_types = TTR.n_types.values

    # fit logarithmic model to TTR curve
    lmodel = LogModel().fit(m_tokens, n_types)
    k_pred = lmodel.predict_k(TTR.m_tokens, normalize=True)
    dim = k_pred.shape[1]
    for n in range(dim):
        TTR[f"lego_{n}"] = TTR[f"lego_{n}"] / TTR.n_types
        TTR[f"lego_{n}_pred"] = k_pred[:, n]

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
    layout = go.Layout(
        title="Type & Legomena Proportions",
        xaxis=dict(title="Tokens"),
        yaxis=dict(title="Proportion of Types", tickformat=",.0%"),
    )
    fig = {"data": data, "layout": layout}
    return fig


# NOTE: to kill "OSError: [Errno 98] Address already in use"
# sudo lsof -t -i tcp:8050 | xargs kill -9

if __name__ == "__main__":
    app.title = "TESTING: " + app.title
    app.run_server(debug=True)
