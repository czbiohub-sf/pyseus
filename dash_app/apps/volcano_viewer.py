import base64
import datetime
from inspect import Attribute
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

from volcano_calculation_layout import calculation_layout
from volcano_plotting_layout import plotting_layout

head, tail = os.path.split(file_dir)
head, tail = os.path.split(head)
sys.path.append(head)

from pyseus import basic_processing as bp
from pyseus import primary_analysis as pa
from pyseus import validation_analysis as va
from pyseus.plotting import plotly_volcano as pv


# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

from dapp import app
from dapp import saved_processed_table, cycle_style_colors, query_panther

# App Layout
layout = html.Div([
    # Header tags
    html.P('Enrichment calculation & Volcano plot generator',
        style={'textAlign': 'center', 'fontSize': 28, 'marginTop': '2%',
            'marginBottom': '1%'}),
    dcc.Tabs(
        id="tabs",
        value='calculation',
        children=[
            dcc.Tab(
                label='Calculate Enrichment & Significance',
                value='calculation',
                children=calculation_layout()
            ),
            dcc.Tab(
                label='Hit-calling & Volcano plot',
                value='plotting',
                children=plotting_layout()
            ),
        ]),
])


@app.callback(
    Output('vol_preloaded_dropdown', 'options'),
    Input('slot_label_1', 'children'),
    Input('slot_label_2', 'children'),
    Input('slot_label_3', 'children'),
    Input('slot_label_4', 'children'),
    Input('slot_label_5', 'children'),
    Input('slot_label_6', 'children'),
)
def load_options(label_1, label_2, label_3, label_4, label_5, label_6):
    """
    automatically load slot labels for the dropdown
    """

    labels = [label_1, label_2, label_3, label_4, label_5, label_6]
    options = []
    for i in np.arange(0, 6):
        option = {'label': 'Slot ' + str(i+1) + ': ' + labels[i], 'value': i}
        options.append(option)
    return options


@app.callback(
    Output('vol_raw_table_upload', 'children'),
    Output('vol_raw_table_upload', 'style'),
    Input('vol_raw_table_upload', 'filename'),
    State('vol_raw_table_upload', 'style'),

    prevent_initial_call=True
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('null_matrix_upload', 'children'),
    Output('null_matrix_upload', 'style'),
    Input('null_matrix_upload', 'filename'),
    State('null_matrix_upload', 'style'),

    prevent_initial_call=True
)
def display_upload_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('samples', 'children'),
    Output('vol_read_table_button', 'style'),
    Output('vol_preload_button', 'style'),
    Input('vol_read_table_button', 'n_clicks'),
    Input('vol_preload_button', 'n_clicks'),
    State('vol_raw_table_upload', 'contents'),


    # preload Input and states

    State('vol_preloaded_dropdown', 'value'),

    State('vol_read_table_button', 'style'),
    State('vol_preload_button', 'style'),

    State('session_id', 'data'),
    prevent_initial_call=True
)
def parse_raw_table(n_clicks, preload_clicks, content, preload_slot,
        button_style, preload_style, session_id):
    """
    group replicates again from the standard file format
    and save the grouped table
    """


    if n_clicks is None and preload_clicks is None:
        raise PreventUpdate

    grouped_table_slot = session_id + 'grouped'
    session_slot = session_id + str(preload_slot)

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]


    if button_id == 'vol_read_table_button':

        # parse txt (csv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0, 1], index_col=0)
        button_style = cycle_style_colors(button_style)

    elif button_id == 'vol_preload_button':

        table = saved_processed_table(session_slot)

        column_tuples = [eval(name) for name in list(table)]
        table.columns = pd.MultiIndex.from_tuples(column_tuples)
        raw_table = table.copy()
        preload_style = cycle_style_colors(preload_style)



    processed_table = raw_table.droplevel(level=0, axis=1).copy()
    processed_table.reset_index(drop=True, inplace=True)
    features = list(raw_table['sample'])
    labels = list(raw_table['metadata'])

    processing = bp.RawTables(file_designated=True)
    processing.transformed_table = processed_table
    processing.sample_cols = features
    processing.info_cols = labels

    processing.group_replicates(reg_exp=r'(.*)_\d+$')
    grouped = processing.grouped_table

    samples = [col[0] for col in list(grouped)]
    samples = list(set(samples))
    samples.remove('metadata')
    samples.sort()

    samples_json = json.dumps(samples)

    _ = saved_processed_table(grouped_table_slot, grouped, overwrite=True)

    return samples_json, button_style, preload_style


@app.callback(
    Output('control_matrix', 'children'),
    Output('null_matrix_button', 'style'),
    Output('control_apply_button', 'style'),
    Input('null_matrix_button', 'n_clicks'),
    Input('control_apply_button', 'n_clicks'),
    Input('samples', 'children'),
    State('null_matrix_upload', 'contents'),

    State('control_matrix', 'children'),
    State('edit_samples', 'value'),
    State('edit_controls', 'value'),

    State('null_matrix_button', 'style'),
    State('control_apply_button', 'style'),

    prevent_intial_call=True
)
def process_control_matrix(upload_clicks, apply_clicks, samples_json,
        matrix_content, matrix, edit_samples, edit_controls, upload_style, apply_style):
    """
    Process or update the manual control sample matrix
    used in significance calculation
    """

    if (samples_json is None) and (upload_clicks is None):
        raise PreventUpdate

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # generate a new null matrix with upload of standard table
    if button_id == 'samples':
        samples = json.loads(samples_json)

        control_matrix = pd.DataFrame()
        # Create a boolean table
        for sample in samples:
            sample_bools = [True if x != sample else False for x in samples]
            control_matrix[sample] = sample_bools

        control_matrix['Samples'] = samples

        control_matrix_json = control_matrix.to_json()

        return control_matrix_json, upload_style, apply_style
        # return control_matrix_json

    # process and update the uploaded custom null matrix
    elif button_id == 'null_matrix_button':

        content_type, content_string = matrix_content.split(',')
        decoded = base64.b64decode(content_string)
        control_matrix = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, index_col=0)

        control_matrix_json = control_matrix.to_json()
        upload_style = cycle_style_colors(upload_style)

        return control_matrix_json, upload_style, apply_style

    # if there are controls to be manually added, apply those changes
    elif button_id == 'control_apply_button':
        matrix = pd.read_json(matrix)
        samples = list(matrix)
        samples.remove('Samples')
        new_controls = [True if sample in edit_controls
            else False for sample in samples]
        for sample in edit_samples:
            matrix[sample] = new_controls

        control_matrix_json = matrix.to_json()
        apply_style = cycle_style_colors(apply_style)

        return control_matrix_json, upload_style, apply_style


@app.callback(
    Output('view_sample_controls', 'options'),
    Output('edit_samples', 'options'),
    Output('edit_controls', 'options'),
    Input('control_matrix', 'children'),
)
def populate_control_options(control_matrix_json):
    if control_matrix_json is None:
        raise PreventUpdate
    else:
        control_matrix = pd.read_json(control_matrix_json)
        samples = list(control_matrix)
        samples.remove('Samples')
        options = []
        for sample in samples:
            option_dict = {'label': sample, 'value': sample}
            options.append(option_dict)

        return options, options, options


@app.callback(
    Output('sample_controls', 'data'),
    Input('view_sample_controls', 'value'),
    Input('control_matrix', 'children'),
    prevent_intial_call=True
)
def review_controls(sample, control_matrix_json):
    if (control_matrix_json is None) or (sample is None):
        raise PreventUpdate
    else:
        control_matrix = pd.read_json(control_matrix_json)

        # Get a list of controls
        controls = control_matrix[['Samples', sample]]
        controls = controls[controls[sample]]
        controls_table = pd.DataFrame(controls[['Samples']])
        controls_table.rename(columns={'Samples': 'sample_control'}, inplace=True)

        return controls_table.to_dict('records')


@app.callback(
    Output('vol_annot_label', 'value'),
    Input('vol_external_annot', 'value'),
    prevent_initial_call=True
)
def fill(annot):
    if annot:
        return annot


@app.callback(
    Output('download_control_matrix', 'data'),
    Output('download_matrix_button', 'style'),
    Input('download_matrix_button', 'n_clicks'),
    State('control_matrix', 'children'),
    State('download_matrix_button', 'style'),
    prevent_initial_call=True
)
def download_matrix(n_clicks, matrix_json, button_style):
    download = pd.read_json(matrix_json)
    download = download.set_index('Samples').reset_index(drop=False)
    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(download.to_csv, 'control_matrix.csv'), button_style


@app.callback(
    Output('calculate_button', 'style'),
    Output('load_button', 'style'),
    Input('calculate_button', 'n_clicks'),
    Input('load_button', 'n_clicks'),
    State('calculation_options', 'value'),
    State('control_matrix', 'children'),
    State('enrichment_option', 'value'),
    State('load_options', 'value'),
    State('significance_table', 'children'),
    State('prep_table_upload', 'contents'),
    State('calculate_button', 'style'),
    State('load_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def calculate_significance(n_clicks, load_clicks, control_opt,
        control_matrix_json, enrichment_opt, load_opt, sig_table_json,
        upload_contents, button_style, load_style, session_id):
    """
    Calculate enrichment and significance based on given options
    or upload enrichment / hits table and cache it, based on user option
    """


    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    grouped_slot = session_id + 'grouped'
    # specify slot for enrichment table
    enriched_slot = session_id + 'enriched'
    hits_slot = session_id + 'hits'
    download_enriched_slot = session_id + 'download'




    # When user uploads the hits table, it is processed by process_hits
    # callback, so do nothing here except change button style

    if (button_id == 'load_button') and (load_opt == 'pre_hits'):
        if upload_contents is not None:

            content_type, content_string = upload_contents.split(',')
            decoded = base64.b64decode(content_string)

            upload_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
                low_memory=False)

            upload_table['fdr'] = upload_table['fdr'].apply(eval)

            _ = saved_processed_table(hits_slot, upload_table, overwrite=True)

            load_style = cycle_style_colors(load_style)

            return button_style, load_style
        else:
            raise PreventUpdate

    # User decides to use calculated enrichment table from calculation tab
    # just need to verify that the enrichment table exists
    elif (button_id == 'load_button') and (load_opt == 'calculated'):

        try:
            # verify that the enrichment table is available
            _ = saved_processed_table(enriched_slot)

        except AttributeError:
            raise PreventUpdate

        load_style = cycle_style_colors(load_style)

        return button_style, load_style

    # if the user is uploading calculated enrichment table, add it to
    # the enrichment table cache
    elif (button_id == 'load_button') and (load_opt == 'pre_enrichment'):
        if upload_contents is None:
            raise PreventUpdate
        content_type, content_string = upload_contents.split(',')
        decoded = base64.b64decode(content_string)

        upload_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False)
        _ = saved_processed_table(enriched_slot, upload_table, overwrite=True)

        load_style = cycle_style_colors(load_style)

        return button_style, load_style

    # user is calculating the enrichment table
    # load grouped table cache and calculate enrichment
    elif button_id == 'calculate_button':
        if enrichment_opt == 'True':
            enrichment_opt = True
        else:
            enrichment_opt = False

        button_style = cycle_style_colors(button_style)
        try:
            # import cached grouped table
            grouped = saved_processed_table(grouped_slot)
        # I need to fix this exception to be specific
        except Exception:
            raise PreventUpdate

        # json reads multi-index tuples as literal strings
        # so convert back to proper multi-level columns
        column_tuples = [eval(name) for name in list(grouped)]
        grouped.columns = pd.MultiIndex.from_tuples(column_tuples)
        grouped.columns = grouped.columns.rename("Samples", level=0)
        grouped.columns = grouped.columns.rename("Replicates", level=1)



        # if control option is manual, import control matrix to assign controls
        if control_opt == 'manual':
            control_matrix = pd.read_json(control_matrix_json)

            analysis = pa.AnalysisTables(imputed_table=grouped,
                exclusion_matrix=control_matrix)
            analysis.simple_pval_enrichment(std_enrich=enrichment_opt)
            analysis.convert_to_standard_table(experiment=False, simple_analysis=True, perseus=False)

            enrichment_table = analysis.standard_hits_table.copy()
            download_table = analysis.simple_pval_table.copy()


        elif control_opt == 'automatic':
            analysis = pa.AnalysisTables(imputed_table=grouped)
            analysis.two_step_bootstrap_pval_enrichment(std_enrich=enrichment_opt)
            analysis.convert_to_standard_table(experiment=False, simple_analysis=False, perseus=False)

            enrichment_table = analysis.standard_hits_table.copy()
            download_table = analysis.two_step_pval_table.copy()

        # process download table for UMAP
        take_cols = []
        for col in list(download_table):
            if col[1] == 'pvals':
                continue
            else:
                take_cols.append(col)

        enrichs = download_table[take_cols].copy()
        new_cols = []
        for col in list(enrichs):
            if col[1] == 'enrichment':
                new_cols.append(('sample', col[0]))
            else:
                new_cols.append(col)

        enrichs.columns = pd.MultiIndex.from_tuples(new_cols)


        _ = saved_processed_table(enriched_slot, enrichment_table, overwrite=True)
        _ = saved_processed_table(download_enriched_slot, enrichs, overwrite=True)

        return button_style, load_style


@app.callback(
    Output('download_umap_table', 'data'),
    Output('download_UMAP_button', 'style'),
    Input('download_UMAP_button', 'n_clicks'),
    State('download_UMAP_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def download_umap(n_clicks, button_style, session_id):
    enriched_slot = session_id + 'download'
    download = saved_processed_table(enriched_slot)
    column_tuples = [eval(name) for name in list(download)]
    download.columns = pd.MultiIndex.from_tuples(column_tuples)
    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(download.to_csv, 'umap_enrichment_table.csv'), button_style


@app.callback(
    Output('download_pval_table', 'data'),
    Output('download_pval_button', 'style'),
    Input('download_pval_button', 'n_clicks'),
    State('download_pval_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def download_pval(n_clicks, button_style, session_id):
    enriched_slot = session_id + 'enriched'
    download = saved_processed_table(enriched_slot)
    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(download.to_csv, 'volcano_enrichment_table.csv'), button_style


@app.callback(
    Output('prep_table_upload', 'children'),
    Output('prep_table_upload', 'style'),
    Input('prep_table_upload', 'filename'),
    State('prep_table_upload', 'style'),

    prevent_initial_call=True
)
def display_upload_name(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('vol_marker_label', 'options'),
    Output('vol_marker_label', 'value'),
    Output('vol_gene_selector', 'options'),
    Output('vol_gene_selector', 'value'),
    Output('vol_merge_key_feature', 'options'),
    Output('vol_merge_key_feature', 'value'),
    Output('enrichment_table_status', 'children'),
    Output('enrichment_table_status', 'style'),
    Input('load_button', 'style'),
    Input('calculate_button', 'style'),
    State('enrichment_table_status', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def check_enrichment_status(load_style, calculate_style, style, session_id):

    enriched_slot = session_id + 'enriched'
    hits_slot = session_id + 'hits'
    enrichment = False
    try:
        sig_table = saved_processed_table(enriched_slot)
        enrichment = True


    except AttributeError:
        try:
            sig_table = saved_processed_table(hits_slot)

        except AttributeError:
            raise PreventUpdate


    meta_cols = [col for col in list(sig_table) if col not in
        ['Unnamed: 0', 'target', 'pvals', 'enrichment']]
    options = []
    label_val = meta_cols[0]
    for col in meta_cols:
        option = {'label': col, 'value': col}
        options.append(option)
        if 'gene' in col.lower():
            label_val = col


    if enrichment:
        style = cycle_style_colors(style)
        return options, label_val, options, label_val, options, label_val,\
            'Enrichment table calculated', style

    else:
        return options, label_val, options, label_val, options, label_val,\
            'Enrichment table not calculated', style


@app.callback(
    Output('hits_button', 'style'),
    Input('hits_button', 'n_clicks'),
    State('prep_table_upload', 'contents'),
    State('thresh_option', 'value'),
    State('fdr', 'value'),
    State('offset', 'value'),
    State('curvature', 'value'),
    State('hits_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True

)
def calculate_hits(hits_clicks, contents, thresh_opt,
        fdr, offset, curvature, hits_style, session_id):
    """
    Calculate the right FDR threshold given options, and call significant hits
    """
    # get the context of the callback trigger

    # set up cache slots
    hits_slot = session_id + 'hits'
    enriched_slot = session_id + 'enriched'

    # load enrichment table from cache
    try:
        sig_table = saved_processed_table(enriched_slot)

    except AttributeError:
        try:
            sig_table = saved_processed_table(hits_slot)

        except AttributeError:
            raise PreventUpdate

    # calculate FDR
    vali = va.Validation(hit_table=sig_table, target_col='target',
        prey_col='prey')
    if thresh_opt == 'hawaii':
        vali.hawaii_fdr(perc=fdr, curvature=curvature, offset_seed=offset,
            experiment=False)
    elif thresh_opt == 'indiv':
        vali.dynamic_fdr(perc=fdr, curvature=curvature, offset_seed=offset,
            experiment=False)

    hits_table = vali.called_table.copy()
    # save hits table to cache
    _ = saved_processed_table(hits_slot, hits_table, overwrite=True)

    hits_style = cycle_style_colors(hits_style)

    return hits_style


@app.callback(
    Output('hits_table_status', 'children'),
    Output('hits_table_status', 'style'),
    Input('hits_button', 'style'),
    Input('load_button', 'style'),
    Input('vol_marker_label', 'value'),
    State('hits_table_status', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def check_hits_status(hits_style, load_style, marker, style, session_id):

    hits_slot = session_id + 'hits'

    try:
        hits_table = saved_processed_table(hits_slot)
    except Exception:
        raise PreventUpdate

    _ = hits_table

    if marker:
        style = cycle_style_colors(style)
        return 'Hits table ready!', style

    else:
        style = cycle_style_colors(style, color_1='#f76868', color_2='#f76868')
        return 'Please select a marker', style


@app.callback(
    Output('download_hits_table', 'data'),
    Output('download_hits_button', 'style'),
    Input('download_hits_button', 'n_clicks'),
    State('download_hits_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def download_hits(n_clicks, button_style, session_id):

    hits_slot = session_id + 'hits'

    download = saved_processed_table(hits_slot)
    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(download.to_csv, 'hits_table.csv'), button_style


@app.callback(
    Output('volcano_dropdown_1', 'options'),
    Input('hits_button', 'style'),
    Input('load_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def populate_volcano_samples(hits_button, loads_button, session_id):
    """
    load sample options in the volcano drop downs
    """

    hits_slot = session_id + 'hits'

    try:
        hits_table = saved_processed_table(hits_slot)
    except Exception:
        raise PreventUpdate

    targets = hits_table['target'].unique()
    targets.sort()
    options = []
    for sample in targets:
        option_dict = {'label': sample, 'value': sample}
        options.append(option_dict)

    return options


@app.callback(
    Output('vol_annot_table_upload', 'children'),
    Output('vol_annot_table_upload', 'style'),
    Input('vol_annot_table_upload', 'filename'),
    State('vol_annot_table_upload', 'style')
)
def display_merge_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(
    Output('vol_merge_key_annot', 'options'),
    Output('vol_external_annot', 'options'),
    Input('vol_annot_table_upload', 'contents'),
    State('vol_annot_table_upload', 'filename')

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
    Output('vol_merge_button', 'style'),
    Input('vol_merge_button', 'n_clicks'),
    State('vol_annot_table_upload', 'contents'),
    State('vol_annot_table_upload', 'filename'),
    State('vol_merge_key_feature', 'value'),
    State('vol_merge_key_annot', 'value'),
    State('vol_external_annot', 'value'),
    State('vol_annot_label', 'value'),
    State('vol_merge_button', 'style'),
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

    session_slot = session_id + 'hits'

    # load cached table
    um_processed_table = saved_processed_table(session_slot)

    # rename keys for proper merge
    annot_table.rename(columns={annot_key: feature_key}, inplace=True)

    merge_table = um_processed_table.merge(annot_table, on=feature_key, how='left')

    drop_subset = list(um_processed_table)
    if 'fdr' in drop_subset:
        drop_subset.remove('fdr')
    merge_table.drop_duplicates(subset=drop_subset, inplace=True)

    rename_label = 'ext_' + annot_label
    merge_table.rename(columns={annot_col: rename_label}, inplace=True)

    _ = saved_processed_table(session_slot, merge_table, overwrite=True)

    return button_style


@app.callback(
    Output('vol_annot_select', 'options'),
    Input('vol_merge_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def fill_ext_options(style, session_id):

    hits_slot = session_id + 'hits'
    try:
        hits_table = saved_processed_table(hits_slot)
    except Exception:
        raise PreventUpdate

    meta_cols = [col for col in list(hits_table) if 'ext_' in col]
    options = []
    for col in meta_cols:
        option = {'label': col, 'value': col}
        options.append(option)

    return options


@app.callback(
    Output('matrix_fig_1', 'figure'),
    Output('estimated_FDR', 'children'),
    Output('preview', 'style'),
    Output('volcano_button_1', 'style'),

    Input('preview', 'n_clicks'),

    Input('volcano_button_1', 'n_clicks'),
    State('plot_options', 'value'),
    State('volcano_dropdown_1', 'value'),
    State('vol_marker_label', 'value'),
    State('vol_annot_select', 'value'),

    State('offset', 'value'),
    State('curvature', 'value'),
    State('estimated_FDR', 'children'),
    State('preview', 'style'),
    State('volcano_button_1', 'style'),



    State('session_id', 'data'),
    prevent_initial_call=True
)
def plot_volcano(volc_click, click_1, checklist, sample, marker, annots,
        offset, curvature, estimate, preview_style, volc_style, session_id):

    hits_slot = session_id + 'hits'
    enriched_slot = session_id + 'enriched'
    enriched = False

    try:
        # import cached grouped table
        hits_table = saved_processed_table(hits_slot)
    except AttributeError:
        try:
            hits_table = saved_processed_table(enriched_slot)
            enriched = True
        except AttributeError:
            raise PreventUpdate

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # creating a preview for seed FDR
    if button_id == 'preview':

        # estimate FDR percentage
        neg_select = hits_table[hits_table['enrichment'] < 0]
        pos_select = hits_table[hits_table['enrichment'] > 0]

        neg_hit = va.hit_count(neg_select, curvature, offset)
        pos_hit = va.hit_count(pos_select, curvature, offset)
        if pos_hit == 0:
            pos_perc = 0
        else:
            pos_perc = np.round(100 * neg_hit / pos_hit, 1)

        estimate = str(pos_perc) + '%'

        # divide the samples to reduce points to plot
        stdev = hits_table['enrichment'].std()
        higher = hits_table[hits_table['enrichment'] > stdev]
        lower = hits_table[hits_table['enrichment'] < -stdev]
        mid = hits_table[~(hits_table.index.isin(higher.index))
            | ~(hits_table.index.isin(lower.index))]

        # sample 5% of the whole data
        mid = mid.sample(frac=0.1, ignore_index=True)
        hits_table = pd.concat([higher, lower, mid])

        # set target name and fdr for Hawaii analysis
        hits_table['target'] = 'hawaii'
        hits_table['interaction'] = False
        hits_table['fdr'] = [[curvature, offset]] * hits_table.shape[0]

        fig = pv.volcano_plot(hits_table, bait='hawaii', marker_mode=False,
            fcd=True, marker=marker, plate='N/A', experiment=False, color=None)

        preview_style = cycle_style_colors(preview_style)

        return fig, estimate, preview_style, volc_style


    # standard-context volcano plotting
    volc_style = cycle_style_colors(volc_style)

    fcd = False
    label = False
    annot = None
    if 'fdr' in checklist:
        fcd = True
        if enriched:
            raise PreventUpdate
    if 'label' in checklist:
        label = True
        if enriched:
            raise PreventUpdate
    if 'ext' in checklist:
        annot = annots

    fig = pv.volcano_plot(hits_table, sample, marker_mode=label, fcd=fcd, marker=marker,
        plate='N/A', experiment=False, color=annot)

    return fig, estimate, preview_style, volc_style


@app.callback(
    Output('vol_selection_count', 'children'),
    Input('matrix_fig_1', 'selectedData'),
    prevent_initial_call=True
)
def print_selection_count(selectedData):
    if selectedData is None:
        PreventUpdate

    num_points = len(selectedData['points'])
    new_str = str(num_points) + ' data points selected'

    return new_str


@app.callback(
    Output('vol_select_button', 'style'),
    Input('vol_select_button', 'n_clicks'),
    State('matrix_fig_1', 'selectedData'),
    State('vol_select_button', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def save_selected_data(n_clicks, selectedData, style, session_id):

    # designate cache ids
    selected_slot = session_id + 'vol_selected'
    complete_slot = session_id + 'hits'

    umap_table = saved_processed_table(complete_slot)

    points = selectedData['points']
    indices = []
    for point in points:
        indices.append(point['customdata'][0])

    selected_table = umap_table[umap_table.index.isin(indices)]
    _ = saved_processed_table(selected_slot, selected_table, overwrite=True)


    style = cycle_style_colors(style)


    return style


@app.callback(
    Output('vol_download_subspace', 'data'),
    Output('vol_download_subspace_button', 'style'),
    Input('vol_download_subspace_button', 'n_clicks'),
    State('vol_download_subspace_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def download_subspace(n_clicks, button_style, session_id):

    umap_slot = session_id + 'vol_selected'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()
        umap_table.reset_index(drop=True, inplace=True)

    except AttributeError:
        raise PreventUpdate

    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(umap_table.to_csv, 'volcano_subspace.csv'), button_style


@app.callback(
    Output('vol_go_top_table', 'data'),
    Output('vol_go_analysis', 'style'),
    Input('vol_go_analysis', 'n_clicks'),
    State('vol_gene_selector', 'value'),
    State('vol_go_cat', 'value'),
    State('vol_pval_cutoff', 'value'),
    State('vol_enrich_cutoff', 'value'),
    State('vol_go_analysis', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def calculate_go(n_clicks, gene_names, category, pval_cutoff, enrch_cutoff, button_style, session_id):

    completed_slot = session_id + 'hits'
    selected_slot = session_id + 'vol_selected'
    go_slot = session_id + 'vol_go'
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
        'pValue': 'pval'
    }, inplace=True)


    return top_ten.to_dict('records'), button_style



@app.callback(
    Output('vol_download_go', 'data'),
    Output('vol_download_go_button', 'style'),
    Input('vol_download_go_button', 'n_clicks'),
    State('vol_download_go_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def download_go(n_clicks, button_style, session_id):

    umap_slot = session_id + 'vol_go'
    if n_clicks is None:
        raise PreventUpdate

    try:
        # verify that the enrichment table is available
        umap_table = saved_processed_table(umap_slot).copy()

    except AttributeError:
        raise PreventUpdate

    button_style = cycle_style_colors(button_style)

    return dcc.send_data_frame(umap_table.to_csv, 'go_analysis.csv'), button_style
