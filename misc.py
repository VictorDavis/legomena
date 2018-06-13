
# bloody dependencies
import math
import plotly.graph_objs as go
import plotly.offline as pyo
import pandas as pd

# Figure 6: tetrakis legomena
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

# Figure 7: fitting different values for k
ttr = pd.read_csv('data/ttr.csv', index_col=1)
df = ttr[ttr['bookid'] == 2701]
tokens = [179402, 90000, 45000, 23000]
for k in [1,2,3]:
    M = tokens[k]
    N = df.loc[M]['types']
    df['types_pred_model{}'.format(k)] = [int(-N*math.log(n/M)*n/M/(1-n/M)) if n != M else N for n in df.index]
data = [
    go.Scatter(
    x = df.index,
    y = df[colname],
    name = displayname,
    mode = 'markers'
) for colname,displayname in {
    'types': 'Types Observed',
    'types_pred_model': 'Types Predicted (M = 180k)',
    'types_pred_model1': 'Types Predicted (M = 90k)',
    'types_pred_model2': 'Types Predicted (M = 45k)',
    'types_pred_model3': 'Types Predicted (M = 23k)'
}.items()]
layout = go.Layout(
    title = 'Type-Token Ratio for k=1,1/2,1/3,1/4',
    xaxis = dict(title = 'Tokens (words)'),
    yaxis = dict(title = 'Types (distinct words)')
)
fig = go.Figure(data = data, layout = layout)
pyo.plot(fig)
