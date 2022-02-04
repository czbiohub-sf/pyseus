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
import os
import sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)

from preprocessing_layout import upload_layout, process_layout

head, tail = os.path.split(file_dir)
head, tail = os.path.split(head)
sys.path.append(head)
from pyseus import basic_processing as bp

from dapp import app
from dapp import saved_processed_table


# App Layout
layout = html.Div([
        # Header tags
        html.P('Data table processing',
            style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'2%',
                'marginBottom': '1%'}),
        dcc.Tabs(
            id="pre_tabs",
            value='raw_table',
            children=[
                dcc.Tab(
                    label='Raw table upload & Feature designation',
                    value='raw_table',
                    children = upload_layout()
                ),
                dcc.Tab(
                    label='Data processing',
                    value='process',
                    children = process_layout()
                ),
                ]),
    ])


@app.callback(
    Output('config_upload', 'children'),
    Output('config_upload', 'style'),
    Input('config_upload', 'filename'),
    State('config_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#DCE7EC'
        return filename, style

@app.callback(
    Output('pp_raw_table_upload', 'children'),
    Output('pp_raw_table_upload', 'style'),
    Input('pp_raw_table_upload', 'filename'),
    State('pp_raw_table_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#DCE7EC'
        return filename, style

@app.callback(
    Output('all_pp_cols', 'children'),
    Output('sample_pp_cols_checklist', 'options'),
    Output('meta_pp_cols_checklist', 'options'),
    Output('read_table_button', 'style'),
    Input('read_table_button', 'n_clicks'),
    State('pp_raw_table_upload', 'contents'),
    State('pp_raw_table_sep', 'value'),
    State('skip_top_rows', 'value'),
    State('read_table_button', 'style'),
    State('session_id', 'data')
    )
def parse_pp_raw_table(n_clicks, content, sep_type, num_skip_row, button_style, session_id):
    """
    initiate QualityControl class with the uploaded proteingroups file
    """

    # create a cache slot for a raw table
    raw_table_slot = session_id + 'raw'
    
    if n_clicks is None:
        raise PreventUpdate
    else:
        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        if num_skip_row == 0:
            skip = None
        else: 
            skip = num_skip_row
        ms_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep_type,
            low_memory=False, skiprows=skip)
        
        # initiate RawTable class
        exp = bp.RawTables(
            experiment_dir='',
            pg_file='',
            intensity_type='Intensity ',
            proteingroup=ms_table,
            file_designated=True)

        saved_processed_table(raw_table_slot, exp.pg_table, overwrite=True)
        all_cols = list(exp.pg_table)
        all_cols_json = json.dumps(all_cols)

        checklist_options = [{'label': col, 'value': col} for col in all_cols]

        if button_style is None:
            button_style = {}
        button_style['background-color'] = '#DCE7EC'

        return all_cols_json, checklist_options,\
            checklist_options, button_style

@app.callback(
    Output('sample_pp_cols_checklist', 'value'),
    Output('apply_sample_re_button', 'style'),    
    Input('apply_sample_re_button', 'n_clicks'),
    State('sample_search_re', 'value'),
    State('all_pp_cols', 'children'),
    State('apply_sample_re_button', 'style'),  
)
def sample_cols_re_search(n_clicks, search_re, all_cols_json, button_style):
    """
    return columns that match regular expression search
    """
    if n_clicks is None:
        raise PreventUpdate
    else:
        button_style['background-color'] = '#DCE7EC'
        all_cols = json.loads(all_cols_json)
        selected_sample_cols = []
        for col in all_cols:
            # check if col name matches with regular expression
            if re.search(search_re, col):
                selected_sample_cols.append(col)
        
        return selected_sample_cols, button_style

@app.callback(
    Output('num_samples_selected', 'children'),
    Input('sample_pp_cols_checklist', 'value')
)
def update_sample_number(checked_list):
    if checked_list is None:
        return ""
    else:
        num_checks = len(checked_list)
        return str(num_checks) + ' sample(s) selected'

@app.callback(
    Output('num_meta_selected', 'children'),
    Input('meta_pp_cols_checklist', 'value')
)
def update_sample_number(checked_list):
    if checked_list is None:
        return ""
    else:
        num_checks = len(checked_list)
        return str(num_checks) + ' sample(s) selected'


@app.callback(
    Output('sample_pp_cols', 'children'),
    Output('meta_pp_cols', 'children'),
    Output('original_sample_pp_cols', 'data'),
    Output('save-cols-button', 'style'),
    Input('save-cols-button', 'n_clicks'),
    State('sample_pp_cols_checklist', 'value'),
    State('meta_pp_cols_checklist', 'value'),
    State('save-cols-button', 'style'),    
)
def select_cols(n_clicks, sample_cols, meta_cols, button_style):
    if n_clicks is None:
        raise PreventUpdate
    else:
        sample_cols_json = json.dumps(sample_cols)
        meta_cols_json = json.dumps(meta_cols)

        col_table = pd.DataFrame()
        col_table['column names'] = sample_cols

        button_style['background-color'] = '#DCE7EC'

        return sample_cols_json, meta_cols_json, col_table.to_dict('records'),\
            button_style

@app.callback(
    Output('renamed_sample_pp_cols', 'data'),
    Output('re_warning', 'children'),
    Output('preview_rename', 'style'),
    Input('preview_rename', 'n_clicks'),
    State('sample_pp_cols', 'children'),
    State('search_re', 'value'),
    State('replacement_re', 'value'),
    State('preview_rename', 'style')
)
def preview_rename_cols(n_clicks, cols_json, search_re, replacement_re, button_style):
    if n_clicks is None:
        raise PreventUpdate

    else:
        search_re = search_re.split(';')
        replacement_re = replacement_re.split(';')
        if len(search_re) != len(replacement_re):
            return None,\
                "Number of search and replacement RE(s) does not match!",\
                {'background-color': '#EF553B'}
        replacement_re = ['' if x == 'NONE' else x for x in replacement_re ]
        intensity_cols = json.loads(cols_json)
        new_cols = bp.sample_rename(intensity_cols,
            search_re, replacement_re)
        
        col_table = pd.DataFrame()
        col_table['new column names'] = new_cols
        button_style['background-color'] = '#DCE7EC'

        return col_table.to_dict('records'), '', button_style

@app.callback(
    Output('process_table_button', 'style'),
    Input('process_table_button', 'n_clicks'),
    State('sample_pp_cols', 'children'),
    State('meta_pp_cols', 'children'),
    State('search_re', 'value'),
    State('replacement_re', 'value'),
    State('filter_rows_check', 'value'),
    State('rename_samples_check', 'value'),
    State('log_transform_check', 'value'),
    State('transform_option', 'value'),
    State('remove_incomplete_rows', 'value'),
    State('merge_reps', 'value'),
    State('merge_option', 'value'),
    State('replace_nulls', 'value'),
    State('replace_option', 'value'),
    State('imputation_dist', 'value'),
    State('imputation_width', 'value'),
    State('tech_reps', 'value'),
    State('process_table_button', 'style'),
    State('session_id', 'data')
)

def process_table(n_clicks, sample_cols_json,
    meta_cols_json, search_re, replacement_re, filter_rows,
    rename_samples, log_transform, transform_opt,
    remove_incomplete_rows, merge_reps, merge_opt, replace_nulls,
    replace_opt, impute_dist, impute_width, tech_reps, button_style, session_id):
    """
    Process table according to user given otions
    """

    # slots for cached tables
    raw_table_slot = session_id + 'raw'
    processed_table_slot = session_id + 'processed'
    
    if n_clicks is None:
        raise PreventUpdate

    # load cached raw table and columns from json formats
    pp_raw_table = saved_processed_table(raw_table_slot)
    sample_cols = json.loads(sample_cols_json)
    meta_cols = json.loads(meta_cols_json)

    # initiate RawTables class from basic_processing
    ms_tables = bp.RawTables(
        experiment_dir=None,
        pg_file='',
        info_cols=meta_cols,
        sample_cols=sample_cols,
        intensity_type='',
        file_designated=True,
        proteingroup=pp_raw_table
    )
    # designate raw table as filtered_table in the class
    # for specific exception to this Dash App
    # and use it as a default output table
    ms_tables.filtered_table = pp_raw_table
    output_table = ms_tables.filtered_table

    # filter rows on MaxQuant contaminants, reverse-seq, and only id by site
    if filter_rows:
        ms_tables.filter_table(select_intensity=False, verbose=False)
        
        output_table = ms_tables.filtered_table

    # rename samples with given search and replacement REs
    if rename_samples:
        search_re = search_re.split(';')
        replacement_re = replacement_re.split(';')
        replacement_re = ['' if x == 'NONE' else x for x in replacement_re ]
        ms_tables.rename_columns(search_re, replacement_re)

        output_table = ms_tables.filtered_table
        
    # transform sample values to log
    if log_transform:
        if transform_opt == 'log2':
            ms_tables.transform_intensities(func=np.log2)
        else:
            ms_tables.transform_intensitites(func=np.log10)
        output_table = ms_tables.transformed_table
    
    # group tables by technical replicates for the following
    # processing options
    if tech_reps == 'yes':
        ms_tables.group_replicates(reg_exp=r'(.*)_\d+$')

    # if options involve any grouping, the output table is defaulted to 
    # grouped table
    if merge_reps or replace_nulls:
        output_table = ms_tables.grouped_table

    # remove rows that do not have at least one sample
    # that has all real values
    if remove_incomplete_rows:
        ms_tables.remove_invalid_rows(verbose=False)
        output_table = ms_tables.preimpute_table
    
    # various null substitution/imputation algorithms
    if replace_nulls:
        if replace_opt == 'zero':
            sample_cols = ms_tables.sample_cols

            # replace null values in sample columns with 0
            output_table[sample_cols] = output_table[sample_cols].fillna(value=0)
        
        elif replace_opt == 'sample_impute':
            ms_tables.bait_impute(distance=impute_dist, width=impute_width)
            output_table = ms_tables.bait_imputed_table
        
        else: 
            ms_tables.prey_impute(distance=impute_dist, width=impute_width)
            output_table = ms_tables.prey_imputed_table
    
    # merge technical replicate with mean or median
    if merge_reps:
        if merge_opt == 'mean':
            mean = True
        else:
            mean=False

        output_table = bp.median_replicates(output_table, mean=mean)
        final_table = bp.dash_output_table(output_table,
            ms_tables.sample_groups, ms_tables.info_cols)

    else: 
        final_table = bp.dash_output_table(output_table,
            ms_tables.sample_cols, ms_tables.info_cols)
    
    # save processed table to cache
    _ = saved_processed_table(processed_table_slot, final_table, overwrite=True)

    # change button color when finished
    button_style['background-color'] = '#DCE7EC'

    return button_style

@app.callback(
    Output('download-ms-table-csv', 'data'),
    Output('download_table', 'style'),
    Input('download_table', 'n_clicks'),
    State('pp_processed_table', 'children'),
    State('ms_save_name', 'value'),
    State('download_table', 'style'),
    State('session_id', 'data'),
    prevent_initial_call=True
)
def download_ms_table(n_clicks, table_json, save_name, button_style, session_id):
    """
    load table from the frontend json format, and save in a csv
    """
    if n_clicks is None:
        raise PreventUpdate

    slot_id = session_id + 'processed'
    download = saved_processed_table(slot_id)

    # json reads multi-index tuples as literal strings
    # so convert back to proper multi-level columns
    column_tuples = [eval(name) for name in list(download)]
    download.columns = pd.MultiIndex.from_tuples(column_tuples)


    button_style['background-color'] = '#DCE7EC'

    return dcc.send_data_frame(download.to_csv,
        save_name), button_style

@app.callback(
    Output('download-configs-csv', 'data'),
    Output('download_configs', 'style'),
    Input('download_configs', 'n_clicks'),
    State('pp_raw_table_sep', 'value'),
    State('skip_top_rows', 'value'),
    State('sample_search_re', 'value'),
    State('search_re', 'value'),
    State('replacement_re', 'value'),
    State('filter_rows_check', 'value'),
    State('rename_samples_check', 'value'),
    State('log_transform_check', 'value'),
    State('transform_option', 'value'),
    State('remove_incomplete_rows', 'value'),
    State('merge_reps', 'value'),
    State('merge_option', 'value'),
    State('replace_nulls', 'value'),
    State('replace_option', 'value'),
    State('imputation_dist', 'value'),
    State('imputation_width', 'value'),
    State('tech_reps', 'value'),
    State('download_configs', 'style'),
    State('configs_save_name', 'value'),
    prevent_initial_call=True
)
def download_configs(n_clicks, table_sep, top_row_skip, sample_search_re,
    search_re, replacement_re, filter_rows,
    rename_samples, log_transform,
    transform_opt, remove_incomplete_rows, merge_reps, merge_opt,
    replace_nulls, replace_opt, impute_dist, impute_width, tech_reps,
    button_style, save_name):

    if n_clicks is None:
        raise PreventUpdate

    configs = pd.DataFrame()
    configs['raw_table_sep'] = [table_sep]
    configs['top_row_skip'] = [top_row_skip]
    configs['sample_search_RE'] = [sample_search_re]
    configs['search_RE'] = [search_re]
    configs['replacement_RE'] = [replacement_re]
    configs['filter_rows'] = [filter_rows]
    configs['rename_samples'] = [rename_samples]
    configs['log_transform'] = [log_transform]
    configs['transform_option'] = [transform_opt]
    configs['remove_incomplete_rows'] = [remove_incomplete_rows]
    configs['merge_replicates'] = [merge_reps]
    configs['merge_option'] = [merge_opt]
    configs['replace_nulls'] = [replace_nulls]
    configs['replace_option'] = [replace_opt]
    configs['impute_distance'] = [impute_dist]
    configs['impute_width'] = [impute_width]
    configs['tech_reps'] = [tech_reps]


    button_style['background-color'] = '#DCE7EC'


    return dcc.send_data_frame(configs.to_csv,
        save_name), button_style

@app.callback(
    Output('pp_raw_table_sep', 'value'),
    Output('skip_top_rows', 'value'),
    Output('sample_search_re', 'value'),
    Output('search_re', 'value'),
    Output('replacement_re', 'value'),
    Output('filter_rows_check', 'value'),
    Output('rename_samples_check', 'value'),
    Output('log_transform_check', 'value'),
    Output('transform_option', 'value'),
    Output('remove_incomplete_rows', 'value'),
    Output('merge_reps', 'value'),
    Output('merge_option', 'value'),
    Output('replace_nulls', 'value'),
    Output('replace_option', 'value'),
    Output('imputation_dist', 'value'),
    Output('imputation_width', 'value'),
    Output('tech_reps', 'value'),
    Output('load-config-button', 'style'),
    Input('load-config-button', 'n_clicks'),
    State('config_upload', 'contents'),
    State('load-config-button', 'style')
)
def preload_configs(n_clicks, content, button_style):
    """automatically update all configs"""
    if n_clicks is None:
        raise PreventUpdate
    else:
        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        
        # values that need to be in a list-form for return
        return_lists = ['filter_rows', 'rename_samples', 'log_transform',
            'remove_incomplete_rows', 'merge_replicates', 'replace_nulls']
        
        # # eval dict for interpreting strings as lists in read_csv
        # converters = {}
        # for col in return_lists:
        #     converters[col] = pd.eval


        # convert uploaded data to csv
        decoded = base64.b64decode(content_string)
        configs = pd.read_csv(io.StringIO(decoded.decode('utf-8')), index_col=0)
        
        for col in return_lists:
            if isinstance(configs[col].item(), str):
                configs[col] = configs[col].apply(pd.eval)
            configs[col] = configs[col].apply(lambda x: x if
                isinstance(x, list) else [])

        button_style = {}
        button_style['background-color'] = '#DCE7EC'
        
        return configs['raw_table_sep'].item(),\
            configs['top_row_skip'].item(),\
            configs['sample_search_RE'].item(),\
            configs['search_RE'].item(),\
            configs['replacement_RE'].item(),\
            configs['filter_rows'].item(),\
            configs['rename_samples'].item(),\
            configs['log_transform'].item(),\
            configs['transform_option'].item(),\
            configs['remove_incomplete_rows'].item(),\
            configs['merge_replicates'].item(),\
            configs['merge_option'].item(),\
            configs['replace_nulls'].item(),\
            configs['replace_option'].item(),\
            configs['impute_distance'].item(),\
            configs['impute_width'].item(),\
            configs['tech_reps'].item(),\
            button_style
 
            
@app.callback(
    Output('pp_processed_table_status', 'data'),
    Input('slot_label_1', 'children'),
    Input('slot_label_2', 'children'),
    Input('slot_label_3', 'children'),
    Input('slot_label_4', 'children'),
    Input('slot_label_5', 'children'),
    Input('slot_label_6', 'children'),
    State('pp_processed_table_status', 'data'),

)
def update_table_status(label_1, label_2, label_3, label_4, label_5, label_6, data):
    table = pd.DataFrame.from_records(data)
    table_names = table['table_name'].to_list()
    table_names[0] = label_1
    table_names[1] = label_2
    table_names[2] = label_3
    table_names[3] = label_4
    table_names[4] = label_5
    table_names[5] = label_6

    table['table_name'] = table_names

    return table.to_dict('records')
        
@app.callback(

    Output('slot_label_4', 'children'),
    Output('slot_label_5', 'children'),
    Output('slot_label_6', 'children'),
    Output('pp_ul_button', 'style'),

    # this is input/states from preprocessing page
    Input('pp_ul_button', 'n_clicks'),
    State('slot_save_name', 'value'),
    State('save_slot', 'value'),


    State('slot_label_4', 'children'),
    State('slot_label_5', 'children'),
    State('slot_label_6', 'children'),
    State('pp_ul_button', 'style'),
    State('session_id', 'data'),

    prevent_initial_call=True
)
def upload_table(n_clicks, pp_save_name, pp_slot,
    label_1, label_2, label_3, button_style, session_id):

    if n_clicks is None:
        raise PreventUpdate

    # since the cache for saved tables start at 4, add 3 to the slot number
    cache_slot = session_id + str(pp_slot + 3)

    # load cached processed table
    processed_slot = session_id + 'processed'

    try: 
        processed_table = saved_processed_table(processed_slot)
    except Exception:
        raise PreventUpdate

    # assignlabels to a list for easy slot designation

    labels = [label_1, label_2, label_3]  
    labels[pp_slot] = pp_save_name

    # save table to designated cache
    _ = saved_processed_table(cache_slot, processed_table, overwrite=True)

    button_style['background-color'] = '#DCE7EC'



    return labels[0], labels[1], labels[2], button_style

if __name__ == "__main__":
    app.run_server(debug=True)