
# bloody dependencies
import plotly.graph_objs as go
import plotly.offline as pyo
import numpy as np
import random

# Experimental observations
ntrials = 1000
deck = range(0, 52)
M_range = list(range(1, 53))
avg_N_suits = np.array([0]*52)
avg_N_valors = np.array([0]*52)
for i in range(1, ntrials+1):
    deal = random.sample(deck, len(deck))
    suits = set([])
    valors = set([])
    N_suits = []
    N_valors = []
    for n in M_range:
        suits.add(deal[n-1] % 4)
        valors.add(deal[n-1] % 13)
        N_suits.append(len(suits))
        N_valors.append(len(valors))
    avg_N_suits = (i-1)*avg_N_suits/i + np.array(N_suits)/i
    avg_N_valors = (i-1)*avg_N_valors/i + np.array(N_valors)/i

# Recursive Formula
Es = [1]
Ev = [1]
for n in range(0, 52):
    Es.append(Es[n]+(52-13*Es[n])/(52-n))
    Ev.append(Ev[n]+(52- 4*Ev[n])/(52-n))

# Continuous Formula
Es2 = []
Ev2 = []
for n in M_range:
    Es2.append(13-13*(1-n/52)**4)
    Ev2.append(4-4*(1-n/52)**13)

data = [
    go.Scatter(
        x = M_range,
        y = avg_N_suits,
        name = "Observed Suits",
        mode = "markers"
    ),
    go.Scatter(
        x = M_range,
        y = avg_N_valors,
        name = "Obvserved Valors",
        mode = "markers"
    ),
    go.Scatter(
        x = M_range,
        y = Es,
        name = "Expected Suits (recursion)",
        mode = "lines"
    ),
    go.Scatter(
        x = M_range,
        y = Ev,
        name = "Expected Valors (recursion)",
        mode = "lines"
    ),
    go.Scatter(
        x = M_range,
        y = Es2,
        name = "Expected Suits (formula)",
        mode = "lines"
    ),
    go.Scatter(
        x = M_range,
        y = Ev2,
        name = "Expected Valors (formula)",
        mode = "lines"
    )
]
layout = go.Layout(
    title = "Playing Deck TTR Curve ({} trials)".format(ntrials),
    xaxis = dict(title = "M (tokens, number of cards dealt)"),
    yaxis = dict(title = "N (types, number of suits/types held)")
)
fig = go.Figure(data = data, layout = layout)
pyo.plot(fig)
