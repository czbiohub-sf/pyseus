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

import numpy as np
import pandas as pd
import plotly.express as px
import re
import sys

from sklearn.cluster import KMeans
import umap

from umap_viewer_layout import create_layout

sys.path.append('../../')
from pyseus import basic_processing as bp
from pyseus.plotting import plotly_umap as pu

# global, immutable variables
transposed_annots = ('sample')

# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# App Layout
app.layout = create_layout()

@app.callback(
    Output('raw_table_upload', 'children'),
    Output('raw_table_upload', 'style'),
    Input('raw_table_upload', 'filename'),
    State('raw_table_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#B6E880'
        return filename, style


@app.callback(
    Output('processed_table', 'children'),
    Output('features', 'children'),
    Output('annots', 'children'),
    Output('table_dims', 'children'), 
    Output('transposed_table', 'children'),
    Output('read_table_button', 'style'),
    Input('read_table_button', 'n_clicks'),
    State('raw_table_upload', 'contents'),
    State('read_table_button', 'style'),
    )
def parse_raw_table(n_clicks, content, button_style):
    """
    initiate QualityControl class with the uploaded proteingroups file
    """

    if n_clicks is None:
        raise PreventUpdate
    else:
        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0,1], index_col=0)
        

        # regular table, features, and annots for UMAP
        processed_table = raw_table.droplevel(level=0, axis=1)
        features = list(raw_table['sample'])
        annots = list(raw_table['metadata'])
        dims = list(raw_table['sample'].shape)

        features_json = json.dumps(features)
        annots_json = json.dumps(annots)
        dims_json = json.dumps(dims)


        # drop duplicates
        processed_table.drop_duplicates(inplace=True)
        processed_table_json = processed_table.to_json()

        # Transposed table, features, and annots for UMAP
        transposed_table = raw_table['sample'].T.reset_index().rename(
            {'index': 'samples'}, axis=1)
        
        transposed_table.drop_duplicates(inplace=True)
        transposed_table_json = transposed_table.to_json()


        if button_style is None:
            button_style = {}
        button_style['background-color'] = '#B6E880'

        return processed_table_json, features_json, annots_json,\
            dims_json, transposed_table_json, button_style

@app.callback(
    Output('feature_dims', 'children'),
    Output('features_checklist', 'options'),
    Output('features_checklist', 'value'),
    Output('label_select', 'options'),
    Output('label_select', 'value'),
    Output('annot_select', 'options'),
    Output('annot_select', 'value'),
    Output('merge_key_feature', 'options'),
    Output('transpose_button', 'style'),
    Input('transpose_button', 'n_clicks'),
    Input('features', 'children'),
    State('annots', 'children'),
    State('table_dims', 'children'),
    State('transpose_button', 'style'),
    prevent_initial_call=True
)
def return_feature_selections(transpose_clicks, features_json, 
    annots_json, dims_json, button_style):

    if transpose_clicks is None:
        transpose_clicks = 0

    # for non transposed option
    if transpose_clicks % 2 == 0:
        features = json.loads(features_json)
        annots = json.loads(annots_json)
        dims = json.loads(dims_json)

        #dimension string formatting
        dim_string = f'{dims[1]} features X {dims[0]} observations, \
            {len(annots)} annotations'
        
        # feature checklist options 
        features_opts = [{'label': feature, 'value': feature}
            for feature in features]

        # labels/annots options
        annot_opts = []
        annot_opts.append({'label': 'None', 'value': 'None'})
        for annot in annots:
            annot_opts.append({'label': annot, 'value': annot})

        annot_val = 'None'
        
        button_style['background-color'] = 'white'
        
        return dim_string, features_opts, features, annot_opts,\
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

        button_style['background-color'] = '#B6E880'

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
        style['background-color'] = '#B6E880'
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
    State('processed_table', 'children'),
    State('transposed_table', 'children'),
    State('merge_key_feature', 'value'),
    State('merge_key_annot', 'value'),
    State('external_annot', 'value'),
    State('merge_button', 'style'),
    prevent_initial_call=True
)
def merge_tables(n_clicks, transpose_clicks, content, filename, processed_table, transposed_table,\
    feature_key, annot_key, annot_label, button_style):

    button_style['background-color'] = '#B6E880'

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
        processed_table = pd.read_json(processed_table)
            
    # for transposed option
    else:
        processed_table = pd.read_json(transposed_table)

    annot_table.rename(columns={annot_key: feature_key}, inplace=True)
    merge_table = processed_table.merge(annot_table, on=feature_key, how='left')
    merge_table.drop_duplicates(subset=list(processed_table), inplace=True)
    
    external_annot = merge_table[annot_label].to_list()
    external_annot_json = json.dumps(external_annot)


    return external_annot_json, button_style


@app.callback(

    Output('umap_fig', 'figure'),
    Output('generate_umap', 'style'),
    Input('generate_umap', 'n_clicks'),
    State('transpose_button', 'n_clicks'),

    State('features_checklist', 'value'),
    State('label_select', 'value'),

    State('annot_options', 'value'),
    State('annot_select', 'value'),
    State('external_annot_series', 'children'),

    State('n_cluster', 'value'),

    State('processed_table', 'children'),
    State('transposed_table', 'children'),

    State('feature_scaling', 'value'),
    State('n_neighbors', 'value'),
    State('min_dist', 'value'),
    State('umap_metric', 'value'),

    State('generate_umap', 'style'),
    
    prevent_initial_call=True
)
def generate_umap(umap_clicks, transpose_clicks, features, label, annot_opts,\
    internal_annot, ext_annot, n_cluster, processed_table, transposed_table, scaling,\
    n_neighbors, min_dist, metric,  button_style):
    """
    Generate umap from all the customizable options
    """

    button_style['background-color'] = '#B6E880'


    if transpose_clicks is None:
        transpose_clicks = 0

    # for non transposed option
    if transpose_clicks % 2 == 0:
        processed_table = pd.read_json(processed_table)
        
    
    # for transposed option
    else:
        processed_table = pd.read_json(transposed_table)
        # manual assignment of features and label
        features = list(processed_table)
        features.remove('samples')
        label = 'samples'
    
    # designate appropriate annotation for colormap
    if annot_opts == 'no_annot':
        annot = 'None'
    elif annot_opts == 'internal':
        annot = internal_annot
    elif annot_opts == 'external':
        external_annots = json.loads(ext_annot)
        processed_table['external'] = external_annots
        annot = 'external'

    processed_table.dropna(subset=features, inplace=True)
    matrix = processed_table[features]

    # scale matrix
    scaled = pu.scale_table(matrix.values, scaling)

    # since clustering annotation must follow dropping Nans, if statement
    # is separately placed here
    if annot_opts == 'cluster':
        clusters = KMeans(n_clusters=n_cluster).fit_predict(scaled)
        clusters = [str(x) for x in clusters]
        processed_table['cluster'] = clusters
        annot = 'cluster'

    # calculate umap dist
    fit = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric
    )
    u = fit.fit_transform(scaled)
    processed_table['umap_1'] = u[: , 0]
    processed_table['umap_2'] = u[:, 1]
    fig = pu.interaction_umap(processed_table, node_name=label, cluster=annot)

    return fig, button_style

    




    


if __name__ == "__main__":
    app.run_server(debug=True)
