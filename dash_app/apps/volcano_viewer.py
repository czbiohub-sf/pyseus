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

from app import app

# App Layout
layout = html.Div([
        # Header tags
        html.P('Enrichment calculation & Volcano plot generator',
            style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'2%',
                'marginBottom': '1%'}),
        dcc.Tabs(
            id="tabs",
            value='calculation',
            children=[
                dcc.Tab(
                    label='Calculate Enrichment & Significance',
                    value='calculation',
                    children = calculation_layout()
                ),
                dcc.Tab(
                    label='Hit-calling & Volcano plot',
                    value='plotting',
                    children = plotting_layout()
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
    
    labels = [label_1, label_2, label_3, label_4, label_5, label_6]
    options = []
    for i in np.arange(0,6):
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
        style['background-color'] = '#DCE7EC'
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
        style['background-color'] = '#DCE7EC'
        return filename, style

@app.callback(
    Output('vol_processed_table', 'children'),
    Output('samples', 'children'),
    Output('vol_read_table_button',   'style'),
    Output('vol_preload_button', 'style'),
    Input('vol_read_table_button', 'n_clicks'),
    Input('vol_preload_button', 'n_clicks'),
    State('vol_raw_table_upload', 'contents'),


    # preload Input and states

    State('vol_preloaded_dropdown', 'value'),
    State('slot_table_1', 'children'),
    State('slot_table_2', 'children'),
    State('slot_table_3', 'children'),
    State('slot_table_4', 'children'),
    State('slot_table_5', 'children'),
    State('slot_table_6', 'children'),

    State('vol_read_table_button', 'style'),
    State('vol_preload_button', 'style'),
    prevent_initial_call=True
    )
def parse_raw_table(n_clicks, preload_clicks, content,  preload_slot,\
        table_1, table_2, table_3, table_4, table_5, table_6, button_style, preload_style):
    """
    group replicates again from the standard file format
    and save the grouped table
    """

  
    if n_clicks is None and preload_clicks is None:
        raise PreventUpdate

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]


    if button_id == 'vol_read_table_button':

        # parse txt (csv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False, header=[0,1], index_col=0)
        if button_style is None:
            button_style = {}
        button_style['background-color'] = '#DCE7EC'
    
    elif button_id == 'vol_preload_button':
        tables = [table_1, table_2, table_3, table_4, table_5, table_6]
        table = pd.read_json(tables[preload_slot])

        column_tuples = [eval(name) for name in list(table)]
        table.columns = pd.MultiIndex.from_tuples(column_tuples)
        raw_table = table.copy()
        if preload_style is None:
            preload_style = {}
        preload_style['background-color'] = '#DCE7EC'


    
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
    
    grouped_json = grouped.to_json()


    return grouped_json, samples_json, button_style, preload_style

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

    if (samples_json == None) and (upload_clicks == None):
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
        upload_style['background-color'] = '#DCE7EC'
        
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
        apply_style['background-color'] = '#DCE7EC'

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
        controls = controls[controls[sample] == True]
        controls_table = pd.DataFrame(controls[['Samples']])
        controls_table.rename(columns={'Samples': 'sample_control'}, inplace=True)
        
        return controls_table.to_dict('records')

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
    button_style['background-color'] = '#DCE7EC'

    return dcc.send_data_frame(download.to_csv, 'control_matrix.csv'), button_style

@app.callback(
    Output('significance_table', 'children'),
    Output('calculate_button', 'style'),
    Output('load_button', 'style'),
    Input('calculate_button', 'n_clicks'),
    Input('load_button', 'n_clicks'),
    State('vol_processed_table', 'children'),
    State('calculation_options', 'value'),
    State('control_matrix', 'children'),
    State('enrichment_option', 'value'),
    State('load_options', 'value'),
    State('significance_table', 'children'),
    State('prep_table_upload', 'contents'),
    State('calculate_button', 'style'),
    State('load_button', 'style'),
    prevent_initial_call=True
)
def calculate_significance(n_clicks, load_clicks, grouped_table_json, control_opt, 
        control_matrix_json, enrichment_opt, load_opt, sig_table_json,
        upload_contents, button_style, load_style):
    """
    Calculate enrichment and significance based on given options
    """

    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if (grouped_table_json is None) and (button_id == 'calculate_button'):
        raise PreventUpdate
    
    elif (button_id == 'load_button') and (load_opt == 'pre_hits'):
        if upload_contents is not None:
            load_style['background-color'] = '#DCE7EC'
            return sig_table_json, button_style, load_style
        else:
            raise PreventUpdate
    
    elif (button_id == 'load_button') and (load_opt == 'calculated'):
       
        if sig_table_json is None:
            raise PreventUpdate

        load_style['background-color'] = '#DCE7EC'
        
        return sig_table_json, button_style, load_style
    
    elif (button_id == 'load_button') and (load_opt == 'pre_enrichment'):
        content_type, content_string = upload_contents.split(',')
        decoded = base64.b64decode(content_string)

        upload_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False)
        upload_table_json = upload_table.to_json()

        load_style['background-color'] = '#DCE7EC'

        return upload_table_json, button_style, load_style

    elif button_id == 'calculate_button':
        if enrichment_opt == 'True':
            enrichment_opt = True
        else:
            enrichment_opt = False

        grouped = pd.read_json(grouped_table_json)
        button_style['background-color'] = '#DCE7EC'

        # json reads multi-index tuples as literal strings
        # so convert back to proper multi-level columns
        column_tuples = [eval(name) for name in list(grouped)]
        grouped.columns = pd.MultiIndex.from_tuples(column_tuples) 
        grouped.columns = grouped.columns.rename("Samples", level=0)
        grouped.columns = grouped.columns.rename("Replicates", level=1) 
        



        if control_opt == 'manual':
            control_matrix = pd.read_json(control_matrix_json)

            analysis = pa.AnalysisTables(imputed_table=grouped,
                exclusion_matrix=control_matrix)
            analysis.simple_pval_enrichment(std_enrich=control_opt)
            analysis.convert_to_standard_table(experiment=False, simple_analysis=True, perseus=False)
            hits_table = analysis.standard_hits_table.copy()
            hits_table_json = hits_table.to_json()          
            
            
            return hits_table_json, button_style, load_style

        elif control_opt =='automatic':
            analysis = pa.AnalysisTables(imputed_table=grouped)
            analysis.two_step_bootstrap_pval_enrichment(std_enrich=enrichment_opt)
            analysis.convert_to_standard_table(experiment=False, simple_analysis=False, perseus=False)

            hits_table = analysis.standard_hits_table.copy()
            hits_table_json = hits_table.to_json()

            return hits_table_json, button_style, load_style

@app.callback(
    Output('download_pval_table', 'data'),
    Output('download_pval_button', 'style'),
    Input('download_pval_button', 'n_clicks'),
    State('significance_table', 'children'),
    State('download_pval_button', 'style'),
    prevent_initial_call=True
)
def download_matrix(n_clicks, table_json, button_style):

    download = pd.read_json(table_json)
    button_style['background-color'] = '#DCE7EC'

    return dcc.send_data_frame(download.to_csv, 'significance_table.csv'), button_style


@app.callback(
    Output('prep_table_upload', 'children'),
    Output('prep_table_upload', 'style'),
    Input('prep_table_upload', 'filename'),
    State('prep_table_upload', 'style'),

    prevent_initial_call=True
)
def display_upload_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#DCE7EC'
        return filename, style


@app.callback(
    Output('vol_marker_label', 'options'),
    Output('enrichment_table_status', 'children'),
    Output('enrichment_table_status', 'style'),
    Input('significance_table', 'children'),
    State('enrichment_table_status', 'style'),

    prevent_initial_call=True
)
def check_enrichment_status(sig_table, style):
    if sig_table is None:
        raise PreventUpdate
    else:
        sig_table = pd.read_json(sig_table)
        
        meta_cols = [col for col in list(sig_table) if col not in ['Unnamed: 0', 'target', 'pvals', 'enrichment']]
        options = []
        for col in meta_cols:
            option = {'label': col, 'value': col}
            options.append(option)

        style['background-color'] = '#DCE7EC'
        return options, 'Enrichment table calculated', style

    
@app.callback(
    Output('hits_table', 'children'),
    Output('hits_button', 'style'),
    Input('hits_button', 'n_clicks'),
    Input('load_button', 'n_clicks'),
    State('load_options', 'value'),
    State('prep_table_upload', 'contents'),
    State('significance_table', 'children'),
    State('thresh_option', 'value'),
    State('fdr', 'value'),
    State('offset', 'value'),
    State('curvature', 'value'),
    State('hits_button', 'style'),

    prevent_initial_call=True

)
def calculate_hits(hits_clicks, load_clicks, load_opt, contents, sig_table_json, thresh_opt,
    fdr, offset, curvature, hits_style):
    """
    Calculate the right FDR threshold given options, and call significant hits
    """
    # get the context of the callback trigger
    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if button_id =='load_button':
        if (contents is None) or (load_opt != 'pre_hits'):
            raise PreventUpdate
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)

        upload_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
            low_memory=False)
        upload_table_json = upload_table.to_json()


        return upload_table_json, hits_style

    elif button_id == 'hits_button':
        if sig_table_json is None:
            raise PreventUpdate
        else:
            sig_table = pd.read_json(sig_table_json)
        vali = va.Validation(hit_table = sig_table, target_col='target', 
            prey_col='prey')
        if thresh_opt == 'hawaii':
            vali.hawaii_fdr(perc=fdr, curvature=curvature, offset_seed=offset,
                experiment=False)
        elif thresh_opt == 'indiv':
            vali.dynamic_fdr(perc=fdr, curvature=curvature, offset_seed=offset,
                experiment=False)
        
        hits_table = vali.called_table.copy()
        hits_table_json = hits_table.to_json()

        hits_style['background-color'] = '#DCE7EC'

        return hits_table_json, hits_style

@app.callback(
    Output('hits_table_status', 'children'),
    Output('hits_table_status', 'style'),
    Input('hits_table', 'children'),
    State('hits_table_status', 'style'),
    prevent_initial_call=True
    )
def check_hits_status(sig_table, style):
    if sig_table is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#DCE7EC'
        return 'Hits table ready!', style

@app.callback(
    Output('download_hits_table', 'data'),
    Output('download_hits_button', 'style'),
    Input('download_hits_button', 'n_clicks'),
    State('hits_table', 'children'),
    State('download_hits_button', 'style'),
    prevent_initial_call=True
)
def download_matrix(n_clicks, table_json, button_style):

    download = pd.read_json(table_json)
    button_style['background-color'] = '#DCE7EC'

    return dcc.send_data_frame(download.to_csv, 'hits_table.csv'), button_style
 
@app.callback(
    Output('volcano_dropdown_1', 'options'),
    Output('volcano_dropdown_2', 'options'),
    Input('hits_table', 'children'),
    prevent_initial_call=True
)
def populate_volcano_samples(hits_table_json):
    if hits_table_json is None:
        raise PreventUpdate
    else:
        hits_table = pd.read_json(hits_table_json)
        targets = hits_table['target'].unique()
        targets.sort()
        options = []
        for sample in targets:
            option_dict = {'label': sample, 'value': sample}
            options.append(option_dict)
        
        return options, options

@app.callback(
    Output('matrix_fig_1', 'figure'),
    Input('volcano_button_1', 'n_clicks'),
    State('hits_table', 'children'),
    State('volcano_dropdown_1', 'value'),
    State('vol_marker_label', 'value'),
    prevent_initial_call=True
)
def plot_volcano(click_1, hits_table_json, sample, marker):

    hits_table = pd.read_json(hits_table_json)
    fig = pv.volcano_plot(hits_table, sample, marker=marker, plate='N/A')

    return fig

@app.callback(
    Output('matrix_fig_2', 'figure'),
    Input('volcano_button_2', 'n_clicks'),
    State('hits_table', 'children'),
    State('volcano_dropdown_2', 'value'),
    State('vol_marker_label', 'value'),
    prevent_initial_call=True)
def plot_volcano(click_1, hits_table_json, sample, marker):

    """
    """
    hits_table = pd.read_json(hits_table_json)
    fig = pv.volcano_plot(hits_table, sample, marker=marker, plate='N/A')

    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
