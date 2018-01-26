"""Functions for generating general quality control plots.

"""
import numpy as np

import plotly.graph_objs as go
from plotly import tools

from genometools.expression import ExpMatrix

from .. import util


def plot_cell_transcript_distribution(
        matrix, name='', color='rgb(31, 119, 180)',
        width=1350, height=600,
        font_size=16, font_family='serif'):
    """Plot the number of transcripts per cell.
    """

    fig = tools.make_subplots(rows=1, cols=2, shared_yaxes=False)
    
    hist_trace = go.Histogram(
        x=matrix.sum(axis=0),
        histnorm='percent',
        marker=dict(
            color=color,
        )
    )

    assert isinstance(matrix, ExpMatrix)
    
    fig.append_trace(hist_trace, 1, 1)
    
    num_transcripts = matrix.sum(axis=0)  # works?
    
    s = np.sort(num_transcripts)[::-1]
    
    step_size = int((matrix.shape[1])/200.0)
    x = s[::step_size]
    y = np.arange(matrix.shape[1])[::step_size] + 1
    
    scatter_trace = go.Scatter(
        x=x,
        y=y,
        mode='lines',
        line=dict(
            color=color,
            width=3.0,
        ),
    )
    
    fig.append_trace(scatter_trace, 1, 2)
    
    layout = go.Layout(
        title = '%s (n=%d)' % (name, matrix.n),
        width=width,
        height=height,
        font=dict(
            size=font_size,
            family=font_family,
        ),
        showlegend=False,
    )
    
    fig['layout'].update(layout)

    xaxis_hist = dict(
        title='Number of transcripts',
    )
    
    yaxis_hist = dict(
        title='Fraction of cells (%)'
    )
    
    fig['layout'].xaxis1.update(xaxis_hist)
    fig['layout'].yaxis1.update(yaxis_hist)
    
    xaxis_scatter = dict(
        title='Transcript threshold',
    )
    
    yaxis_scatter = dict(
        title=('Number of cells above threshold'),
        autorange=False,
        range=[0, matrix.shape[1]*1.05],
    )
    
    fig['layout'].xaxis2.update(xaxis_scatter)
    fig['layout'].yaxis2.update(yaxis_scatter)
    
    return fig


def plot_transcriptome_components(
        matrix, species='human', name='',
        width=950, height=800, font_size=16, font_family='serif'):
    """Plots showing mitochondrial and ribosomal transcriptome components.
    
    TODO: docstring
    """

    def generate_scatter_trace(matrix, sel_genes, color):

        transcripts = matrix.sum(axis=0)
        sel_sum = matrix.loc[sel_genes].sum(axis=0)
        sel_frac = sel_sum / matrix.sum(axis=0)

        trace = go.Scatter(
            x=transcripts,
            y=100*sel_frac,
            text=matrix.cells,
            mode='markers',
            marker=dict(
                opacity=0.7,
                color=color,
            ),
        )

        return trace

    def generate_hist_trace(matrix, sel_genes, color):
        
        transcripts = matrix.sum(axis=0)
        sel_sum = matrix.loc[sel_genes].sum(axis=0)
        sel_frac = sel_sum / matrix.sum(axis=0)

        trace = go.Histogram(
            y=100*sel_frac,
            autobiny=False,
            ybins=dict(
                start=0,
                end=100.1,
                size=5.0001,
            ),
            marker=dict(
                color=color,
            ),
            histnorm='percent',
        )

        return trace
    
    fig = tools.make_subplots(rows=2, cols=2, shared_yaxes=True)

    mito_genes = util.get_mitochondrial_genes(species=species)
    ribo_genes = util.get_ribosomal_genes(species=species)

    mito_color = 'rgb(255, 127, 14)'
    ribo_color = 'rgb(31, 119, 180)'

    try:
        mito_trace1 = generate_hist_trace(matrix, mito_genes, mito_color)
        mito_trace2 = generate_scatter_trace(
            matrix, mito_genes, mito_color)
    except KeyError:
        pass
    else:
        fig.append_trace(mito_trace1, 1, 1)
        fig.append_trace(mito_trace2, 1, 2)

    try:
        ribo_trace1 = generate_scatter_trace(
            matrix, ribo_genes, ribo_color)
        ribo_trace2 = generate_hist_trace(matrix, ribo_genes, ribo_color)
    except KeyError:
        pass
    else:
        fig.append_trace(ribo_trace1, 2, 1)
        fig.append_trace(ribo_trace2, 2, 2)

    if name:
        name = name + ' '
    title = '%s(n=%d)' % (name, matrix.shape[1])

    layout = go.Layout(
        title=title,
        width=width,
        height=height,
        font=dict(
            size=font_size,
            family='font_family',
        ),
        showlegend=False,
        bargap=0.1,
    )
    #fig['layout'].update(height=600, width=600, title='i <3 subplots')

    xaxis_hist = dict(
        title='Fraction of cells (%)'
    )

    xaxis_scatter = dict(
        title='Number of transcripts'
    )

    yaxis_mito = dict(
        title='Fraction of mitochondrial<br> transcripts (%)',
        autorange=False,
        range=[0, 100],
        zeroline=True,
    )

    yaxis_ribo = dict(
        title='Fraction of ribosomal<br> transcripts (%)',
        autorange=False,
        range=[0, 100],
        zeroline=True,
    )

    fig['layout'].update(layout)

    fig['layout'].xaxis1.update(xaxis_hist)
    fig['layout']['xaxis1']['title'] = ''
    fig['layout'].xaxis3.update(xaxis_hist)

    #fig['layout'].xaxis.update(xaxis)
    fig['layout'].xaxis2.update(xaxis_scatter)
    fig['layout']['xaxis2']['title'] = ''
    fig['layout'].xaxis4.update(xaxis_scatter)

    fig['layout'].yaxis1.update(yaxis_mito)
    fig['layout'].yaxis2.update(yaxis_ribo)

    #fig['layout']['yaxis2'].update(yaxis_mito)
    #fig['layout']['yaxis2']['title'] = ''

    #fig['layout']['yaxis4'].update(yaxis_ribo)
    #fig['layout']['yaxis4']['title'] = ''

    return fig
