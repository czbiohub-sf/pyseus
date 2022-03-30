import base64
import datetime
import markdown
import io
import json
import simplejson
import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dcc
from dash import html
from dash import dash_table

from plotly import graph_objs as go

import numpy as np
import pandas as pd
import plotly.express as px
import re
import os
import sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
from sklearn.cluster import KMeans
import umap

from matrix_viewer_layout import calculation_layout, plotting_layout

head, tail = os.path.split(file_dir)
head, tail = os.path.split(head)
sys.path.append(head)
from pyseus import basic_processing as bp
from pyseus.plotting import plotly_umap as pu
from pyseus.plotting import plotly_heatmap as ph

from dapp import app
from dapp import saved_processed_table, cycle_style_colors

# global, immutable variables
transposed_annots = ('sample')

# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']





# App Layout
layout = html.Div([
    # Header tags
    html.P('Clustergram generator',
        style={'textAlign': 'center', 'fontSize': 28, 'marginTop': '2%',
            'marginBottom': '1%'}),
    dcc.Tabs(
        id="tabs",
        value='options',
        children=[
            dcc.Tab(
                label='Clustergram options',
                value='options',
                children=calculation_layout()
            ),
            dcc.Tab(
                label='Plot Clustergram',
                value='plotting',
                children=plotting_layout()
            ),
        ]),
])


@app.callback(
    Output('mat_preloaded_dropdown', 'options'),
    Input('slot_label_1', 'children'),
    Input('slot_label_2', 'children'),
    Input('slot_label_3', 'children'),
    Input('slot_label_4', 'children'),
    Input('slot_label_5', 'children'),
    Input('slot_label_6', 'children'),
)
def load_options(label_1, label_2, label_3, label_4, label_5, label_6):
    """
    automatically populate slot labels based on the upload/save list
    """
    labels = [label_1, label_2, label_3, label_4, label_5, label_6]
    options = []
    for i in np.arange(0, 6):
        option = {'label': 'Slot ' + str(i+1) + ': ' + labels[i], 'value': i}
        options.append(option)
    return options


@app.callback(
    Output('mat_raw_table_upload', 'children'),
    Output('mat_raw_table_upload', 'style'),
    Input('mat_raw_table_upload', 'filename'),
    State('mat_raw_table_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('features_checklist', 'options'),
    Output('features_checklist', 'value'),
    Output('label_select', 'options'),
    Output('index_select', 'options'),
    Output('data_metrics', 'data'),
    Output('color_button', 'n_clicks'),
    Output('mat_read_table_button', 'style'),
    Output('mat_preload_button', 'style'),

    # upload Input and states
    Input('mat_read_table_button', 'n_clicks'),
    Input('mat_preload_button', 'n_clicks'),
    State('mat_raw_table_upload', 'contents'),
    State('mat_read_table_button', 'style'),
    State('color_button', 'n_clicks'),

    # preload Input and states

    State('mat_preloaded_dropdown', 'value'),
    State('mat_preload_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def parse_raw_table(n_clicks, preload_clicks, content, button_style, color_clicks,
        preload_slot, preload_style, session_id):

    """
    Load cached table or upload a new table, and cache it to specific
    clustergram-designated slot.
    """

    if n_clicks is None and preload_clicks is None:
        raise PreventUpdate

    # combine unique session ID with designated slot id for loading cache
    session_slot = session_id + str(preload_slot)

    # unique cache id for the processed table used for clust page
    clust_slot = session_id + 'clust'


    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if button_id == 'mat_read_table_button':
        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0, 1], index_col=0)

        button_style = cycle_style_colors(button_style)

    elif button_id == 'mat_preload_button':

        table = saved_processed_table(session_slot)

        column_tuples = [eval(name) for name in list(table)]
        table.columns = pd.MultiIndex.from_tuples(column_tuples)
        raw_table = table.copy()
        preload_style = cycle_style_colors(preload_style)


    # regular table, features, and annots for matrix
    processed_table = raw_table.droplevel(level=0, axis=1).copy()


    features = list(raw_table['sample'])
    features.sort()
    labels = list(raw_table['metadata'])
    labels.sort()

    # feature checklist options
    features_opts = [{'label': feature, 'value': feature}
        for feature in features]


    # labels/annots options
    label_opts = []
    label_opts.append({'label': 'None', 'value': 'None'})
    for label in labels:
        label_opts.append({'label': label, 'value': label})

    _ = saved_processed_table(clust_slot, processed_table, overwrite=True)



    if color_clicks is None:
        color_clicks = 1
    else:
        color_clicks += 1

    matrix = processed_table[features]
    metrics = [{
        'min': [np.round(matrix.values.min(), 2)],
        'max': [np.round(matrix.values.max(), 2)],
        'avg': [np.round(matrix.values.mean(), 2)],
        'stdev': [np.round(matrix.values.std(), 2)]
    }]

    return features_opts, features,\
        label_opts, label_opts, metrics, color_clicks, button_style, preload_style



@app.callback(
    Output('color_bar', 'figure'),
    Input('scale_data_button', 'n_clicks'),
    Input('color_button', 'n_clicks'),
    State('colorscale_min', 'value'),
    State('colorscale_max', 'value'),
    State('colormap', 'value')
)
def generate_colormap(scale_data_clicks, color_clicks,
        min, max, colormap):

    fig, _ = ph.color_map(min, max, colors=colormap)
    return fig


@app.callback(
    Output('matrix_fig', 'figure'),
    Output('generate_matrix', 'style'),

    Input('generate_matrix', 'n_clicks'),
    State('features_checklist', 'value'),
    State('label_select', 'value'),
    State('index_select', 'value'),

    State('colorscale_min', 'value'),
    State('colorscale_max', 'value'),
    State('colormap', 'value'),

    State('cluster_checks', 'value'),
    State('tick_checks', 'value'),

    State('generate_matrix', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def generate_clustergram(n_clicks, features, label, index,
        zmin, zmax, colormap, cluster_checks, tick_checks, button_style, session_id):
    """
    returns plotly figure of cluster heatmap
    """

    if n_clicks is None:
        raise PreventUpdate

    button_style = cycle_style_colors(button_style)

    # read cache of cluster table
    clust_id = session_id + 'clust'
    processed_table = saved_processed_table(clust_id)


    # generate the color map
    _, hexmap = ph.color_map(zmin, zmax, colormap)


    # default bait clustering variables
    bait_leaves = None
    bait_clust = False

    # cluster samples
    if 'bait_clust' in cluster_checks:
        bait_leaves = ph.bait_leaves(processed_table, features, grouped=False,
            verbose=False)
        bait_clust = True

    prey_leaves = ph.prey_leaves(processed_table, features, index_id=index, grouped=False,
        verbose=False)

    heatmap = ph.dendro_heatmap(processed_table, prey_leaves, hexmap,
        zmin, zmax, label, features, index_id=index, bait_leaves=bait_leaves, bait_clust=bait_clust,
        verbose=False)

    x_tick = False
    y_tick = False
    if 'sample_ticks' in tick_checks:
        x_tick = True
    if 'obs_ticks' in tick_checks:
        y_tick = True

    layout = go.Layout(
        xaxis={'showticklabels': x_tick},
        yaxis={'showticklabels': y_tick})

    fig = go.Figure(data=heatmap, layout=layout)

    return fig, button_style


if __name__ == "__main__":
    app.run_server(debug=True)
