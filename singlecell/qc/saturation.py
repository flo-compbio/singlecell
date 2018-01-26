import math

import numpy as np
import plotly.graph_objs as go


def plot_saturation(matrix, fractions=None, bin_edges=None):
    """Plots the sequencing saturation for different fractions of reads.
    
    TODO: docstring"""

    #assert isinstance(matrix, SparseExpMatrix)

    if fractions is None:
        fractions = np.arange(0.05, 0.99, 0.05)

    #X = matrix.values
    p, n = matrix.shape

    X = matrix.values

    # calculate total number of transcripts per sample
    num_transcripts = (X>0).sum(axis=0, dtype=np.uint32)

    #num_transcripts = np.sum(S.array() > 0, axis=0, dtype=np.uint32)
    
    S = np.ones((fractions.size, n), dtype=np.float64)
    
    for i, frac in enumerate(fractions):
        assert 0 <= frac <= 1.0
        
        if frac == 1.0:
            continue
        
        # calculate inclusion probability for each transcript (UMI)
        # => this is the probability that at least one read of this transcript
        #    is included in the analysis
        
        # this unnecessarily creates a copy - very wasteful use of memory

        # calculate inclusion probability for all transcripts of each sample
        sample_sat = np.zeros(n, dtype=np.float64)
        for j in range(n):
            v = X[:, j]
            if not isinstance(v, np.ndarray):
                # for compatibility with scipy.sparse matrices
                v = v.toarray()
            #print(v[:10])
            sample_sat[j] = np.sum(1.0 - np.power(1.0-frac, v))
        
        sample_sat = sample_sat / np.float64(num_transcripts)
        S[i, :] = sample_sat

    #print(S[:3, :10])

    data = []

    #print(num_transcripts[:5])

    bin_auto = False
    if bin_edges is None:
        bin_auto = True
        bin_edges = np.power(10, np.arange(7))
    d = np.digitize(num_transcripts, bins=bin_edges)
    #bc = np.bincount(d)
    #for b, num_cells in enumerate(bc):
    for b in range(1, min(np.amax(d)+1, len(bin_edges))):
        sel = (d == b).nonzero()[0]
        if sel.size > 0:
            if bin_auto:
                name = '10<sup>%d</sup> - 10<sup>%d</sup> (n=%d)' \
                    % (b-1, b, sel.size)
            else:
                name = '%d - %d (n=%d)' \
                       % (bin_edges[b-1], bin_edges[b], sel.size)

            mean_sat = np.mean(S[:, sel], axis=1)

            trace = go.Scatter(
                x=100*fractions,
                y=100*mean_sat,
                mode='lines+markers',
                line=dict(
                    width=5,
                ),
                marker=dict(
                    size=10,
                ),
                name=name,
            )

            data.append(trace)

    layout = go.Layout(
        #title = 'Total transcripts: %d' % exp_matrix.iloc[:, 0].sum(),
        margin=dict(
            l=100,
        ),
        xaxis=dict(
            title='Subsampling rate (%)',
        ),
        yaxis=dict(
            title='Fraction of transcripts<br>recovered (%)',
            autorange=False,
            range=[0, 100],
        ),
        font=dict(
            size=20,
            family='Serif',
        ),
        showlegend=True,
    )

    fig = go.Figure(data=data, layout=layout)
    return fig
