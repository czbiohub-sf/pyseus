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

from flask_caching.backends import FileSystemCache

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

from umap_viewer_layout import plotting_layout, customize_layout

from flask_caching import Cache


head, tail = os.path.split(file_dir)
head, tail = os.path.split(head)
sys.path.append(head)
from pyseus import basic_processing as bp
from pyseus.plotting import plotly_umap as pu

# global, immutable variables
transposed_annots = ('sample')

from dapp import app
from dapp import saved_processed_table


# App Layout
layout = html.Div([
        # Header tags
        html.P('The UMAP generator',
            style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'2%',
                'marginBottom': '1%'}),
        dcc.Tabs(
            id="tabs",
            value='calculation',
            children=[
                dcc.Tab(
                    label='Customize features and annotations',
                    value='calculation',
                    children = customize_layout()
                ),
                dcc.Tab(
                    label='Plot UMAP',
                    value='plotting',
                    children = plotting_layout()
                ),
                ],
            style={'marginBottom':'2%'}),
    ])

@app.callback(
    Output('um_preloaded_dropdown', 'options'),
    Input('slot_label_1', 'children'),
    Input('slot_label_2', 'children'),
    Input('slot_label_3', 'children'),
    Input('slot_label_4', 'children'),
    Input('slot_label_5', 'children'),
    Input('slot_label_6', 'children'),
)
def load_options(label_1, label_2, label_3, label_4, label_5, label_6):
    
    labels = [label_1, label_2, label_3, label_4, label_5, label_6]
    options = []
    for i in np.arange(0,6):
        option = {'label': 'Slot ' + str(i+1) + ': ' + labels[i], 'value': i}
        options.append(option)
    return options 


@app.callback(
    Output('um_raw_table_upload', 'children'),
    Output('um_raw_table_upload', 'style'),
    Input('um_raw_table_upload', 'filename'),
    State('um_raw_table_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#DCE7EC'
        return filename, style
    

@app.callback(
    Output('um_features', 'children'),
    Output('um_annots', 'children'),
    Output('table_dims', 'children'), 
    Output('um_read_table_button', 'style'),
    Output('um_preload_button', 'style'),
    Input('um_read_table_button', 'n_clicks'),
    Input('um_preload_button', 'n_clicks'),

    State('um_raw_table_upload', 'contents'),

    # preload Input and states

    State('um_preloaded_dropdown', 'value'),
    State('um_read_table_button', 'style'),
    State('um_preload_button', 'style'),
    State('session_id', 'data')

    )
def parse_um_raw_table(n_clicks, preload_clicks, content, 
    preload_slot, button_style, preload_style, session_id):
    """
    initiate QualityControl class with the uploaded proteingroups file
    """

    if n_clicks is None and preload_clicks is None:
        raise PreventUpdate

    session_slot = session_id + str(preload_slot)
    um_slot = session_id + 'umap'

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if button_id == 'um_read_table_button':

        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        um_raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0,1], index_col=0)
        if button_style is None:
            button_style = {}
        button_style['background-color'] = '#DCE7EC'        

    elif button_id == 'um_preload_button':
        table = saved_processed_table(session_slot)

        column_tuples = [eval(name) for name in list(table)]
        table.columns = pd.MultiIndex.from_tuples(column_tuples)
        um_raw_table = table.copy()
        if preload_style is None:
            preload_style = {}
        preload_style['background-color'] = '#DCE7EC'

    # regular table, features, and annots for UMAP
    um_processed_table = um_raw_table.droplevel(level=0, axis=1)
    um_features = list(um_raw_table['sample'])
    annots = list(um_raw_table['metadata'])
    dims = list(um_raw_table['sample'].shape)

    um_features_json = json.dumps(um_features)
    annots_json = json.dumps(annots)
    dims_json = json.dumps(dims)


    # drop duplicates
    um_processed_table.drop_duplicates(inplace=True)
    saved_processed_table(um_slot, um_processed_table)

    # Transposed table, features, and annots for UMAP
    transposed_table = um_raw_table['sample'].T.reset_index().rename(
        {'index': 'samples'}, axis=1)
    
    transposed_table.drop_duplicates(inplace=True)
    trans_slot = session_id + 'transposed'
    saved_processed_table(trans_slot, transposed_table)

    return um_features_json, annots_json,\
        dims_json, button_style, preload_style

@app.callback(
    Output('feature_dims', 'children'),
    Output('um_features_checklist', 'options'),
    Output('um_features_checklist', 'value'),
    Output('um_label_select', 'options'),
    Output('um_label_select', 'value'),
    Output('annot_select', 'options'),
    Output('annot_select', 'value'),
    Output('merge_key_feature', 'options'),
    Output('transpose_button', 'style'),
    Input('transpose_button', 'n_clicks'),
    Input('um_features', 'children'),
    State('um_annots', 'children'),
    State('table_dims', 'children'),
    State('transpose_button', 'style'),
    prevent_initial_call=True
)
def return_feature_selections(transpose_clicks, um_features_json, 
    annots_json, dims_json, button_style):

    if um_features_json is None and transpose_clicks is None:
        raise PreventUpdate

    if transpose_clicks is None:
        transpose_clicks = 0

    # for non transposed option
    if transpose_clicks % 2 == 0:
        um_features = json.loads(um_features_json)
        annots = json.loads(annots_json)
        dims = json.loads(dims_json)

        #dimension string formatting
        dim_string = f'{dims[1]} features X {dims[0]} observations, \
            {len(annots)} annotations'
        
        # feature checklist options 
        um_features_opts = [{'label': feature, 'value': feature}
            for feature in um_features]

        # labels/annots options
        annot_opts = []
        annot_opts.append({'label': 'None', 'value': 'None'})
        for annot in annots:
            annot_opts.append({'label': annot, 'value': annot})

        annot_val = 'None'
        
        button_style['background-color'] = 'white'
        
        return dim_string, um_features_opts, um_features, annot_opts,\
            annot_val, annot_opts, annot_val, annot_opts, button_style
    
    # for transposed option
    else:
        dims = json.loads(dims_json)
        dim_string = f'{dims[0]} features X {dims[1]} observations, \
            0 annotations'

        feature_str = 'Feature selection N/A'
        feature_opts = [{'label': feature_str, 'value': 'None'}]
        feature_val = ['None']

        label_opts = [{'label': 'sample names', 'value': 'sample'}]

        annot_str = 'No internal annotation available\
            on transposed table'
        annot_opts = [{'label': annot_str, 'value': 'None'}]
        annot_val = 'None'

        button_style['background-color'] = '#DCE7EC'

        return dim_string, feature_opts, feature_val, label_opts, annot_val,\
            annot_opts, annot_val, label_opts, button_style


@app.callback(
    Output('annot_table_upload', 'children'),
    Output('annot_table_upload', 'style'),
    Input('annot_table_upload', 'filename'),
    State('annot_table_upload', 'style')
)
def display_merge_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#DCE7EC'
        return filename, style

@app.callback(
    Output('merge_key_annot', 'options'),
    Output('external_annot', 'options'),
    Input('annot_table_upload', 'contents'),
    State('annot_table_upload', 'filename')

)
def fill_external_keys(content, filename):
    if content is None:
        raise PreventUpdate
    else:
        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        if 'csv' in filename:
            annot_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif 'tsv' in filename:
            annot_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
        
        cols = list(annot_table)
        # feature checklist options 
        opts = [{'label': feature, 'value': feature}
            for feature in cols if "Unnamed" not in feature]
        
        return opts, opts  



@app.callback(
    Output('external_annot_series', 'children'),
    Output('merge_button', 'style'),
    Input('merge_button', 'n_clicks'),
    State('transpose_button', 'n_clicks'),
    State('annot_table_upload', 'contents'),
    State('annot_table_upload', 'filename'),
    State('merge_key_feature', 'value'),
    State('merge_key_annot', 'value'),
    State('external_annot', 'value'),
    State('merge_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def merge_tables(n_clicks, transpose_clicks, content, filename, \
    feature_key, annot_key, annot_label, button_style, session_id):

    if n_clicks is None:
        raise PreventUpdate

    button_style['background-color'] = '#DCE7EC'

    # parse txt (tsv) file as pd df from upload
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)

    if 'csv' in filename:
        annot_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    elif 'tsv' in filename:
        annot_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
        

    if transpose_clicks is None:
        transpose_clicks = 0

    # for non transposed option
    if transpose_clicks % 2 == 0:
        session_slot = session_id + 'umap'

    # for transposed option
    else:
        session_slot = session_id + 'transposed'
    
    um_processed_table = saved_processed_table(session_slot)

    annot_table.rename(columns={annot_key: feature_key}, inplace=True)
    merge_table = um_processed_table.merge(annot_table, on=feature_key, how='left')
    merge_table.drop_duplicates(subset=list(um_processed_table), inplace=True)
    
    external_annot = merge_table[annot_label].to_list()
    external_annot_json = json.dumps(external_annot)


    return external_annot_json, button_style


@app.callback(

    Output('umap_fig', 'figure'),
    Output('generate_umap', 'style'),
    Input('generate_umap', 'n_clicks'),
    State('transpose_button', 'n_clicks'),

    State('um_features_checklist', 'value'),
    State('um_label_select', 'value'),

    State('annot_options', 'value'),
    State('annot_select', 'value'),
    State('external_annot_series', 'children'),

    State('n_cluster', 'value'),


    State('feature_scaling', 'value'),
    State('n_neighbors', 'value'),
    State('min_dist', 'value'),
    State('umap_metric', 'value'),

    State('generate_umap', 'style'),
    State('session_id', 'data'),
    
    prevent_initial_call=True
)
def generate_umap(umap_clicks, transpose_clicks, um_features, label, annot_opts,\
    internal_annot, ext_annot, n_cluster, scaling,\
    n_neighbors, min_dist, metric,  button_style, session_id):
    """
    Generate umap from all the customizable options
    """
    if umap_clicks is None:
        raise PreventUpdate

    button_style['background-color'] = '#DCE7EC'


    if transpose_clicks is None:
        transpose_clicks = 0

    # for non transposed option
    if transpose_clicks % 2 == 0:
        session_slot = session_id + 'umap'
        um_processed_table = saved_processed_table(session_slot)

    # for transposed option
    else:
        session_slot = session_id + 'transposed'
        um_processed_table = saved_processed_table(session_slot)


        # manual assignment of features and label
        um_features = list(um_processed_table)
        um_features.remove('samples')
        label = 'samples'
    
    # designate appropriate annotation for colormap
    if annot_opts == 'no_annot':
        annot = 'None'
    elif annot_opts == 'internal':
        annot = internal_annot
    elif annot_opts == 'external':
        external_annots = json.loads(ext_annot)
        um_processed_table['external'] = external_annots
        annot = 'external'

    um_processed_table.dropna(subset=um_features, inplace=True)
    matrix = um_processed_table[um_features]

    # scale matrix
    scaled = pu.scale_table(matrix.values, scaling)

    # since clustering annotation must follow dropping Nans, if statement
    # is separately placed here
    if annot_opts == 'cluster':
        clusters = KMeans(n_clusters=n_cluster).fit_predict(scaled)
        clusters = [str(x) for x in clusters]
        um_processed_table['cluster'] = clusters
        annot = 'cluster'

    # calculate umap dist
    fit = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric
    )
    u = fit.fit_transform(scaled)
    um_processed_table['umap_1'] = u[: , 0]
    um_processed_table['umap_2'] = u[:, 1]
    fig = pu.interaction_umap(um_processed_table, node_name=label, cluster=annot)

    complete_slot = session_id +'completed'
    um_processed_table = saved_processed_table(complete_slot, um_processed_table)

    return fig, button_style

    
@app.callback(
    Output('download_umap', 'data'),
    Output('download_umap_button', 'style'),
    Input('download_umap_button', 'n_clicks'),
    State('download_umap_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def download_matrix(n_clicks, button_style, session_id):
    if n_clicks is None:
        raise PreventUpdate
    slot = session_id + 'completed'
    download = saved_processed_table(slot)
    button_style['background-color'] = '#DCE7EC'

    return dcc.send_data_frame(download.to_csv, 'umap.csv'), 

    


if __name__ == "__main__":
    app.run_server(debug=True)
