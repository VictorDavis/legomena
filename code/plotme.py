
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.offline as pyo
import plotly.graph_objs as go

# all have to be read into memory?!
books = pd.read_csv('data/books.csv', index_col = 0)
ttr   = pd.read_csv('data/ttr.csv')

# Fig Data: Real vs. Optimum Type/Token scatterplot
def graphData():
    data = [
        go.Scatter(
            x = books['tokens'],
            y = books['types'],
            name = 'Actual',
            text = books['title'],
            mode = 'markers'
        ),
        go.Scatter(
            x = books['Mz'],
            y = books['Nz'],
            name = 'Optimum',
            text = books['title'],
            mode = 'markers'
        )
    ]
    layout = go.Layout(
        title = 'Actual vs Optimum Type/Token Values',
        xaxis = dict(title = 'Corpus Size (Tokens)'),
        yaxis = dict(title = 'Vocabulary Size (Types)')
    )
    fig = {'data': data, 'layout': layout}
    return fig

# Fig RMSE: RMSE by corpus size
def graphRMSE():
    data = [
        go.Scatter(
            x = books['tokens'],
            y = books[col],
            name = label,
            text = books['title'],
            mode = 'markers'
        )
        for col,label in {'RMSE_heaps':'Heap\'s Law','RMSE_taylor':'Taylor Series','RMSE_model':'Logarithmic Model'}.items()
    ]
    layout = go.Layout(
        title = 'Root Mean Square Error (% of Types)',
        xaxis = dict(title = 'Corpus Size (Tokens)'),
        yaxis = dict(title = 'RMSE %')
    )
    fig = {'data': data, 'layout': layout}
    return fig

# Show plots in dashboard
app = dash.Dash()

app.layout = html.Div(children = [
    dcc.Graph(id = 'figData', figure = graphData()), # actual vs. optimum types by corpus size
    dcc.Graph(id = 'figRMSE', figure = graphRMSE()), # RMSE by corpus size
    dcc.Dropdown(
        id = 'book-selector',
        options = [
            { 'label': '{}: {}'.format(books.loc[key]['author'], books.loc[key]['title']), 'value': key } for key in books.index
        ],
        value = 2701 # default to Moby Dick
    ),
    dcc.Graph(id = 'figTTR'), # TTR Curve
    dcc.Graph(id = 'figLego') # Hapax & n-legomena proportions
])

# Fig TTR: Actual vs Predicted types = f(tokens)
@app.callback(
    Output('figTTR', 'figure'),
    [ Input('book-selector', 'value') ]
)
def graph1(bookid):
    df = ttr[ttr['bookid'] == bookid]
    title = books.loc[bookid]['title']
    Mz    = books.loc[bookid]['Mz']
    data = [
        go.Scatter(
        x = df['tokens'],
        y = df[colname],
        name = displayname,
        mode = 'lines' if '_pred' in colname else 'markers'
    ) for colname,displayname in {
        'types': 'Types Observed',
        'types_pred_heaps': 'Types Predicted (Heap\'s Law)',
        'types_pred_taylor': 'Types Predicted (Taylor Series)',
        'types_pred_model': 'Types Predicted (Logarithmic Model)',
        'hapax': 'Hapaxes Observed',
        'hapax_pred_model': 'Hapaxes Predicted',
        'dis': 'Dis Legomena Observed',
        'dis_pred_model': 'Dis Predicted',
        'tris': 'Tris Legomena Observed',
        'tris_pred_model': 'Tris Predicted',
        'tetrakis': 'Tetrakis Legomena Observed',
        'tetrakis_pred_model': 'Tetrakis Predicted',
        'pentakis': 'Pentakis Legomena Observed',
        'pentakis_pred_model': 'Pentakis Predicted'
    }.items()]
    layout = go.Layout(
        title = 'Type-Token Ratio for {} (Mz={})'.format(title, Mz),
        xaxis = dict(title = 'Tokens'),
        yaxis = dict(title = 'Types')
    )
    fig = {'data': data, 'layout': layout}
    return fig

# Fig Lego: Hapax & n-legomena proportions as sample size grows
@app.callback(
    Output('figLego', 'figure'),
    [ Input('book-selector', 'value') ]
)
def graphLego(bookid):
    df = ttr[ttr['bookid'] == bookid]
    title = books.loc[bookid]['title']
    data = [go.Scatter(
        x = df['tokens'],
        y = df[colname],
        name = displayname,
        mode = 'lines' if '_pred' in colname else 'markers'
    ) for colname,displayname in {
        'hapax_frac': 'Hapaxes Observed',
        'hapax_frac_pred': 'Hapaxes Predicted',
        'dis_frac': 'Dis Legomena Observed',
        'dis_frac_pred': 'Dis Predicted',
        'tris_frac': 'Tris Legomena Observed',
        'tris_frac_pred': 'Tris Predicted',
        'tetrakis_frac': 'Tetrakis Legomena Observed',
        'tetrakis_frac_pred': 'Tetrakis Predicted',
        'pentakis_frac': 'Pentakis Legomena Observed',
        'pentakis_frac_pred': 'Pentakis Predicted'
    }.items()]
    layout = go.Layout(
        title = 'Legomena Proportions for {}'.format(title),
        xaxis = dict(title = 'Tokens'),
        yaxis = dict(title = 'Fraction of Types')
    )
    fig = {'data': data, 'layout': layout}
    return fig

if __name__ == '__main__':
    app.run_server()
