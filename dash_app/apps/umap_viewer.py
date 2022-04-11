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
import time


head, tail = os.path.split(file_dir)
head, tail = os.path.split(head)
sys.path.append(head)
from pyseus import basic_processing as bp
from pyseus.plotting import plotly_umap as pu

# global, immutable variables
transposed_annots = ('sample')

from dapp import app
from dapp import saved_processed_table, cycle_style_colors, query_panther, collapsible_style


# App Layout
layout = html.Div([
    # Header tags
    html.P('The UMAP generator',
        style={'textAlign': 'center', 'fontSize': 28, 'marginTop': '2%',
            'marginBottom': '1%'}),
    dcc.Tabs(
        id="tabs",
        value='calculation',
        children=[
            dcc.Tab(
                label='Customize features and annotations',
                value='calculation',
                children=customize_layout()
            ),
            dcc.Tab(
                label='Plot UMAP',
                value='plotting',
                children=plotting_layout()
            ),
        ],
        style={'marginBottom': '2%'}),
])


@app.callback(
    Output('label_div', 'style'),
    Output('label_section', 'children'),
    Input('label_section', 'n_clicks'),
    State('label_section', 'children'),
    State('label_div', 'style'),
    prevent_initial_call=True
)
def process_collapse(n_clicks, button_txt, section_1):
    sections = [section_1]
    button_txt, sections = collapsible_style(button_txt, sections)

    return sections[0], button_txt


@app.callback(
    Output('go_div', 'style'),
    Output('go_section', 'children'),
    Input('go_section', 'n_clicks'),
    State('go_section', 'children'),
    State('go_div', 'style'),
    prevent_initial_call=True
)
def go_collapse(n_clicks, button_txt, section_1):
    sections = [section_1]
    button_txt, sections = collapsible_style(button_txt, sections)

    return sections[0], button_txt


@app.callback(
    Output('feature_div', 'style'),
    Output('feature_section', 'children'),
    Input('feature_section', 'n_clicks'),
    State('feature_section', 'children'),
    State('feature_div', 'style'),
    prevent_initial_call=True
)
def feature_collapse(n_clicks, button_txt, section_1):
    sections = [section_1]
    button_txt, sections = collapsible_style(button_txt, sections)

    return sections[0], button_txt


@app.callback(
    Output('umops_div', 'style'),
    Output('umops_section', 'children'),
    Input('umops_section', 'n_clicks'),
    State('umops_section', 'children'),
    State('umops_div', 'style'),
    prevent_initial_call=True
)
def umops_collapse(n_clicks, button_txt, section_1):
    sections = [section_1]
    button_txt, sections = collapsible_style(button_txt, sections)

    return sections[0], button_txt


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
    Output('um_raw_table_upload', 'children'),
    Output('um_raw_table_upload', 'style'),
    Input('um_raw_table_upload', 'filename'),
    State('um_raw_table_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('umap_table_upload', 'children'),
    Output('umap_table_upload', 'style'),
    Input('umap_table_upload', 'filename'),
    State('umap_table_upload', 'style')
)
def display_upload_umap_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
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
    Load cached table or upload a new table, and cache it to specific
    UMAP-designated slot.
    """

    if n_clicks is None and preload_clicks is None:
        raise PreventUpdate

    # combine unique session ID with designated slot id for loading cache
    session_slot = session_id + str(preload_slot)

    # unique cache id for the processed table used for UMAP page
    um_slot = session_id + 'umap'

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Condition for processing uploaded table
    if button_id == 'um_read_table_button':

        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        um_raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0, 1], index_col=0)
        button_style = cycle_style_colors(button_style)

    # retrieving table from cache in specific slot
    elif button_id == 'um_preload_button':
        table = saved_processed_table(session_slot)

        # json processing for multi indexing
        column_tuples = [eval(name) for name in list(table)]
        table.columns = pd.MultiIndex.from_tuples(column_tuples)
        um_raw_table = table.copy()

        preload_style = cycle_style_colors(preload_style)

    # processing table, features, and annots for UMAP
    um_processed_table = um_raw_table.droplevel(level=0, axis=1)
    um_features = list(um_raw_table['sample'])
    um_features.sort()
    annots = list(um_raw_table['metadata'])
    annots.sort()

    dims = list(um_raw_table['sample'].shape)

    # save the data as JSON for client-side storage
    um_features_json = json.dumps(um_features)
    annots_json = json.dumps(annots)
    dims_json = json.dumps(dims)


    # drop duplicates and save table to cache
    um_processed_table.drop_duplicates(inplace=True)
    saved_processed_table(um_slot, um_processed_table, overwrite=True)

    # Transposed table, features, and annots for UMAP
    transposed_table = um_raw_table['sample'].T.reset_index().rename(
        {'index': 'samples'}, axis=1)

    transposed_table.drop_duplicates(inplace=True)
    trans_slot = session_id + 'transposed'
    # save to transpose-specific cache
    saved_processed_table(trans_slot, transposed_table, overwrite=True)

    return um_features_json, annots_json,\
        dims_json, button_style, preload_style


@app.callback(
    Output('feature_dims', 'children'),
    Output('um_features_checklist', 'options'),
    Output('um_features_checklist', 'value'),
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
    """
    From the data extracted from processed table, populate options for
    features, labels, annots, etc.
    """

    if um_features_json is None and transpose_clicks is None:
        raise PreventUpdate

    if transpose_clicks is None:
        transpose_clicks = 0

    # for non transposed option
    if transpose_clicks % 2 == 0:
        um_features = json.loads(um_features_json)
        annots = json.loads(annots_json)
        dims = json.loads(dims_json)

        # dimension string formatting
        dim_string = f'{dims[1]} features X {dims[0]} observations, \
            {len(annots)} annotations'

        # feature checklist options
        um_features_opts = [{'label': feature, 'value': feature}
            for feature in um_features]


        button_style['background-color'] = 'white'

        return dim_string, um_features_opts, um_features, button_style

    # for transposed option
    else:
        dims = json.loads(dims_json)
        dim_string = f'{dims[0]} features X {dims[1]} observations, \
            0 annotations'

        feature_str = 'Feature selection N/A'
        feature_opts = [{'label': feature_str, 'value': 'None'}]
        feature_val = ['None']


        button_style = cycle_style_colors(button_style)

        return dim_string, feature_opts, feature_val, button_style


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
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('merge_key_annot', 'options'),
    Output('external_annot', 'options'),
    Input('annot_table_upload', 'contents'),
    State('annot_table_upload', 'filename')

)
def fill_external_keys(content, filename):
    """
    populate dropdown options from uploaded annotation table
    """
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
    Output('merge_button', 'style'),
    Input('merge_button', 'n_clicks'),
    State('annot_table_upload', 'contents'),
    State('annot_table_upload', 'filename'),
    State('merge_key_feature', 'value'),
    State('merge_key_annot', 'value'),
    State('external_annot', 'value'),
    State('annot_label', 'value'),
    State('merge_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def merge_tables(n_clicks, content, filename,
        feature_key, annot_key, annot_col, annot_label, button_style, session_id):
    """
    Load umap table from cache, and merge it with external annotation table.
    Save merged series to client-side.
    """

    if n_clicks is None:
        raise PreventUpdate

    button_style = cycle_style_colors(button_style)

    # parse txt (tsv) file as pd df from upload
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)

    if 'csv' in filename:
        annot_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    elif 'tsv' in filename or 'txt' in filename:
        annot_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')

    annot_table = annot_table[[annot_key, annot_col]]

    session_slot = session_id + 'completed'

    # load cached table
    um_processed_table = saved_processed_table(session_slot)

    # rename keys for proper merge
    annot_table.rename(columns={annot_key: feature_key}, inplace=True)

    merge_table = um_processed_table.merge(annot_table, on=feature_key, how='left')

    merge_table.drop_duplicates(subset=list(um_processed_table), inplace=True)

    rename_label = 'ext_' + annot_label
    merge_table.rename(columns={annot_col: rename_label}, inplace=True)

    _ = saved_processed_table(session_slot, merge_table, overwrite=True)

    return button_style


@app.callback(
    Output('x_select', 'options'),
    Output('y_select', 'options'),
    Input('final_features', 'children'),
    prevent_initial_call=True
)
def populate_data_opts(final_features_json):

    features = json.loads(final_features_json)
    default_opts = [
        {'label': 'UMAP 1', 'value': 'umap_1'},
        {'label': 'UMAP 2', 'value': 'umap_2'}]

    for feature in features:
        default_opts.append({'label': feature, 'value': feature})

    return default_opts, default_opts


@app.callback(
    Output('cluster_button', 'style'),
    Input('cluster_button', 'n_clicks'),
    State('n_cluster', 'value'),
    State('final_features', 'children'),
    State('cluster_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def make_clusters(n_clicks, num_clust, um_features_json, button_style, session_id):

    session_slot = session_id + 'completed'

    button_style = cycle_style_colors(button_style)

    # load cached table
    um_processed_table = saved_processed_table(session_slot)

    um_features = json.loads(um_features_json)

    matrix = um_processed_table[um_features]

    # scale matrix
    scaled = pu.scale_table(matrix.values, 'standard')
    # K means
    clusters = KMeans(n_clusters=num_clust).fit_predict(scaled)
    clusters = [str(x) for x in clusters]


    um_processed_table['cluster_' + str(num_clust)] = clusters

    _ = saved_processed_table(session_slot, um_processed_table, overwrite=True)

    return button_style


@app.callback(
    Output('generate_umap', 'style'),
    Output('umap_load_button', 'style'),
    Output('final_features', 'children'),
    Input('generate_umap', 'n_clicks'),

    Input('umap_load_button', 'n_clicks'),
    State('umap_table_upload', 'contents'),
    State('umap_load_button', 'style'),

    State('transpose_button', 'n_clicks'),
    State('um_annots', 'children'),
    State('um_features_checklist', 'value'),
    State('feature_scaling', 'value'),
    State('n_neighbors', 'value'),
    State('min_dist', 'value'),
    State('umap_metric', 'value'),
    State('random_state', 'value'),

    State('generate_umap', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def generate_umap(umap_clicks, upload_clicks, content, upload_style, transpose_clicks, um_annots_json,
        um_features, scaling, n_neighbors, min_dist, metric, random_state,
        generate_style, session_id):
    """
    Generate umap from all the customizable options
    """
    if umap_clicks is None and upload_clicks is None:
        raise PreventUpdate



    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Condition for processing uploaded table
    if button_id == 'umap_load_button':

        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        um_raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0, 1], index_col=0)

        upload_style = cycle_style_colors(upload_style)

        um_features = list(um_raw_table['sample'])
        um_features.sort()
        um_features_json = json.dumps(um_features)

        umap_table = um_raw_table.droplevel(level=0, axis=1)
        session_slot = session_id + 'completed'

        _ = saved_processed_table(session_slot, umap_table, overwrite=True)

        return generate_style, upload_style, um_features_json


    if transpose_clicks is None:
        transpose_clicks = 0

    # cached load for non transposed option
    if transpose_clicks % 2 == 0:
        session_slot = session_id + 'umap'
        um_processed_table = saved_processed_table(session_slot)

        annots = json.loads(um_annots_json)
        um_processed_table = um_processed_table[um_features + annots]

    # cached load for transposed option
    else:
        session_slot = session_id + 'transposed'
        um_processed_table = saved_processed_table(session_slot)


        # manual assignment of features and label
        um_features = list(um_processed_table)
        um_features.remove('samples')


    _ = list(um_processed_table)
    um_processed_table.dropna(subset=um_features, inplace=True)
    matrix = um_processed_table[um_features]


    # scale matrix
    scaled = pu.scale_table(matrix.values, scaling)


    # configure random state
    if random_state == 'None':
        # calculate umap dist
        fit = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric
        )
    else:
        random_state = int(random_state)
        fit = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state
        )
    u = fit.fit_transform(scaled)
    um_processed_table['umap_1'] = u[: , 0]
    um_processed_table['umap_2'] = u[:, 1]



    # save table with umap coordinates for able cache
    complete_slot = session_id + 'completed'
    um_processed_table = saved_processed_table(complete_slot, um_processed_table, overwrite=True)

    um_features_json = json.dumps(um_features)
    generate_style = cycle_style_colors(generate_style)



    return generate_style, upload_style, um_features_json


@app.callback(
    Output('gene_selector', 'options'),
    Output('gene_selector', 'value'),
    Output('merge_key_feature', 'options'),
    Output('merge_key_feature', 'value'),
    Output('um_label_select', 'options'),
    Output('um_label_select', 'value'),
    Output('annot_select', 'options'),
    Output('annot_select', 'value'),
    Input('generate_umap', 'style'),
    Input('umap_load_button', 'style'),
    Input('merge_button', 'style'),
    Input('cluster_button', 'style'),
    State('final_features', 'children'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def populate_options(input_1, input_2, input_3, input_4, features_json, session_id):

    umap_slot = session_id + 'completed'

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    features = json.loads(features_json)

    cols = list(umap_table)
    annots = [col for col in cols if col not in features and 'umap_' not in col]
    labels = [col for col in annots if 'cluster_' not in col and 'ext_' not in col]


    # labels/annots options
    annot_opts = []
    annot_opts.append({'label': 'None', 'value': 'None'})
    for annot in annots:
        annot_opts.append({'label': annot, 'value': annot})

    annot_val = 'None'

    label_val = 'None'
    # labels/annots options
    label_opts = []
    label_opts.append({'label': 'None', 'value': 'None'})
    for label in labels:
        label_opts.append({'label': label, 'value': label})
        if 'gene' in label.lower():
            label_val = label

    return label_opts, label_val, label_opts, label_val, label_opts, label_val, annot_opts, annot_val


@app.callback(
    Output('umap_fig', 'figure'),
    Output('plot_button', 'style'),
    Output('search_button', 'style'),
    Input('plot_button', 'n_clicks'),
    Input('search_button', 'n_clicks'),
    State('search_plot', 'value'),
    State('search_button', 'style'),
    State('um_label_select', 'value'),
    State('annot_select', 'value'),
    State('marker_color', 'value'),
    State('opacity', 'value'),
    State('x_select', 'value'),
    State('y_select', 'value'),
    State('plot_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def plot_umap(n_clicks, search_clicks, search_term, search_style, label, annot,
        marker_color, opacity, x_val, y_val, button_style, session_id):

    umap_slot = session_id + 'completed'
    if n_clicks is None and search_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if button_id == 'search_button':
        search_lower = search_term.lower()
        umap_table[search_term] = umap_table[label].map(lambda x: search_term if search_lower
            in str(x).lower() else None)

        # umap generation
        fig = pu.interaction_umap(umap_table, node_name=label, cluster=search_term,
            unlabelled_color=marker_color, unlabelled_opacity=opacity * 0.5, x=x_val, y=y_val,
            unlabelled_hover=False, search=True)
        search_style = cycle_style_colors(search_style)


    else:
        # umap generation
        fig = pu.interaction_umap(umap_table, node_name=label, cluster=annot,
            unlabelled_color=marker_color, unlabelled_opacity=opacity, x=x_val, y=y_val)

        button_style = cycle_style_colors(button_style)

    if 'umap' not in x_val:
        fig.add_vline(x=0, line_width=1)
        fig.add_hline(y=0, line_width=1)




    return fig, button_style, search_style


@app.callback(
    Output('annot_label', 'value'),
    Input('external_annot', 'value'),
    prevent_initial_call=True
)
def fill(annot):
    if annot:
        return annot


@app.callback(
    Output('download_umap', 'data'),
    Output('download_umap_button', 'style'),
    Input('download_umap_button', 'n_clicks'),
    State('um_features', 'children'),
    State('final_features', 'children'),
    State('download_umap_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def download_umap(n_clicks, orig_features_json, final_features_json, button_style, session_id):

    umap_slot = session_id + 'completed'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    if orig_features_json:
        features = json.loads(orig_features_json)
    elif final_features_json:
        features = json.loads(final_features_json)

    column_tups = []
    for col in list(umap_table):
        if col in features:
            column_tups.append(('sample', col))
        else:
            column_tups.append(('metadata', col))

    umap_table.columns = pd.MultiIndex.from_tuples(column_tups)

    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(umap_table.to_csv, 'umap_table.csv'), button_style


@app.callback(
    Output('annot_table_dl', 'data'),
    Output('annot_dl_button', 'style'),
    Input('annot_dl_button', 'n_clicks'),
    State('um_features', 'children'),
    State('final_features', 'children'),
    State('annot_dl_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def download_annot_umap(n_clicks, orig_features_json, final_features_json, button_style, session_id):

    umap_slot = session_id + 'completed'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    if orig_features_json:
        features = json.loads(orig_features_json)
    elif final_features_json:
        features = json.loads(final_features_json)

    column_tups = []
    for col in list(umap_table):
        if col in features:
            column_tups.append(('sample', col))
        else:
            column_tups.append(('metadata', col))

    umap_table.columns = pd.MultiIndex.from_tuples(column_tups)

    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(umap_table.to_csv, 'annotated_umap_table.csv'), button_style


@app.callback(
    Output('umap_table_status', 'children'),
    Output('umap_table_status', 'style'),
    Input('generate_umap', 'style'),
    Input('umap_load_button', 'style'),
    State('umap_table_status', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def check_umap_status(hits_style, load_style, style, session_id):

    hits_slot = session_id + 'completed'

    try:
        hits_table = saved_processed_table(hits_slot)
    except Exception:
        raise PreventUpdate

    _ = hits_table

    style = cycle_style_colors(style)
    return 'UMAP table ready!', style


@app.callback(
    Output('selection_count', 'children'),
    Output('selected_data', 'data'),
    Input('umap_fig', 'selectedData'),
    State('um_label_select', 'value'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def print_selection_count(selectedData, label, session_id):
    if selectedData is None:
        PreventUpdate

    num_points = len(selectedData['points'])
    new_str = str(num_points) + ' data points selected'


    # designate cache ids
    selected_slot = session_id + 'selected'
    complete_slot = session_id + 'completed'

    umap_table = saved_processed_table(complete_slot)

    points = selectedData['points']
    indices = []
    for point in points:
        indices.append(point['customdata'][0])

    selected_table = umap_table[umap_table.index.isin(indices)]
    selected_table.reset_index(drop=True, inplace=True)

    labels = selected_table.rename(columns={label: 'marker'})[['marker']].sort_values(
        by='marker')


    _ = saved_processed_table(selected_slot, selected_table, overwrite=True)

    return new_str, labels.to_dict('records')


@app.callback(
    Output('download_subspace', 'data'),
    Output('download_subspace_button', 'style'),
    Input('download_subspace_button', 'n_clicks'),
    State('download_subspace_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def download_subspace(n_clicks, button_style, session_id):

    umap_slot = session_id + 'selected'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(umap_table.to_csv, 'umap_subspace.csv'), button_style


@app.callback(
    Output('go_top_table', 'data'),
    Output('go_analysis', 'style'),
    Input('go_analysis', 'n_clicks'),
    State('gene_selector', 'value'),
    State('go_cat', 'value'),
    State('pval_cutoff', 'value'),
    State('enrich_cutoff', 'value'),
    State('go_analysis', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def calculate_go(n_clicks, gene_names, category, pval_cutoff, enrch_cutoff, button_style, session_id):

    completed_slot = session_id + 'completed'
    selected_slot = session_id + 'selected'
    go_slot = session_id + 'umap_go'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        completed = saved_processed_table(completed_slot).copy()
        selected = saved_processed_table(selected_slot).copy()

    except AttributeError:
        raise PreventUpdate

    all_genes = completed[gene_names].apply(
        lambda x: str(x).upper().split(';')).explode().drop_duplicates().to_list()
    selected_genes = selected[gene_names].apply(
        lambda x: str(x).upper().split(';')).explode().drop_duplicates().to_list()



    go_table = query_panther(selected_genes, all_genes, pval_thresh=pval_cutoff,
        enrichment_thresh=enrch_cutoff, biological=category)

    _ = saved_processed_table(go_slot, go_table, overwrite=True)

    button_style = cycle_style_colors(button_style)

    top_ten = go_table.iloc[:10]

    top_ten['expected'] = top_ten['expected'].apply(lambda x: np.round(x, 2))
    top_ten['fold_enrichment'] = top_ten['fold_enrichment'].apply(lambda x: np.round(x, 2))

    top_ten.rename(columns={
        'go_term_label': 'go_term',
        'number_in_list': 'num',
        'fold_enrichment': 'fold',
        'pValue': 'pval'}, inplace=True)


    return top_ten.to_dict('records'), button_style


@app.callback(
    Output('download_go', 'data'),
    Output('download_go_button', 'style'),
    Input('download_go_button', 'n_clicks'),
    State('download_go_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def download_go(n_clicks, button_style, session_id):

    umap_slot = session_id + 'umap_go'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(umap_table.to_csv, 'umap_go_analysis.csv'), button_style
