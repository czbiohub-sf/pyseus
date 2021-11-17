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

from preprocessing_layout import create_layout

sys.path.append('../../')
from pyseus import basic_processing as bp


# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# App Layout
app.layout = create_layout()


@app.callback(
    Output('raw_table_filename', 'children'),
    Output('raw_table_filename', 'style'),
    Input('raw_table_upload', 'filename'),
    State( 'raw_table_filename', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#B6E880'
        return filename, style


@app.callback(
    Output('config_filename', 'children'),
    Output('config_filename', 'style'),
    Input('config_upload', 'filename'),
    State('config_filename', 'style'),
)
def display_upload_config_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style['background-color'] = '#B6E880'
        return filename, style

@app.callback(
    Output('raw_table', 'children'),
    Output('all_cols', 'children'),
    Output('sample_cols_checklist', 'options'),
    Output('meta_cols_checklist', 'options'),
    Output('read_table_button', 'style'),
    Input('read_table_button', 'n_clicks'),
    State('raw_table_upload', 'contents'),
    State('raw_table_sep', 'value'),
    State('skip_top_rows', 'value'),
    State('read_table_button', 'style'),
    )
def parse_raw_table(n_clicks, content, sep_type, num_skip_row, button_style):
    """
    initiate QualityControl class with the uploaded proteingroups file
    """

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

        filtered_table = exp.pg_table.to_json()
        all_cols = list(exp.pg_table)
        all_cols_json = json.dumps(all_cols)

        checklist_options = [{'label': col, 'value': col} for col in all_cols]

        if button_style is None:
            button_style = {}
        button_style['background-color'] = '#B6E880'

        return filtered_table, all_cols_json, checklist_options,\
            checklist_options, button_style

@app.callback(
    Output('sample_cols_checklist', 'value'),
    Output('apply_sample_re_button', 'style'),    
    Input('apply_sample_re_button', 'n_clicks'),
    State('sample_search_re', 'value'),
    State('all_cols', 'children'),
    State('apply_sample_re_button', 'style'),  
)
def sample_cols_re_search(n_clicks, search_re, all_cols_json, button_style):
    """
    return columns that match regular expression search
    """
    if n_clicks is None:
        raise PreventUpdate
    else:
        button_style['background-color'] = '#B6E880'
        all_cols = json.loads(all_cols_json)
        selected_sample_cols = []
        for col in all_cols:
            # check if col name matches with regular expression
            if re.search(search_re, col):
                selected_sample_cols.append(col)
        
        return selected_sample_cols, button_style

@app.callback(
    Output('num_samples_selected', 'children'),
    Input('sample_cols_checklist', 'value')
)
def update_sample_number(checked_list):
    if checked_list is None:
        return ""
    else:
        num_checks = len(checked_list)
        return str(num_checks) + ' sample(s) selected'

@app.callback(
    Output('num_meta_selected', 'children'),
    Input('meta_cols_checklist', 'value')
)
def update_sample_number(checked_list):
    if checked_list is None:
        return ""
    else:
        num_checks = len(checked_list)
        return str(num_checks) + ' sample(s) selected'


@app.callback(
    Output('sample_cols', 'children'),
    Output('meta_cols', 'children'),
    Output('original_sample_cols', 'data'),
    Output('save-cols-button', 'style'),
    Input('save-cols-button', 'n_clicks'),
    State('sample_cols_checklist', 'value'),
    State('meta_cols_checklist', 'value'),
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

        button_style['background-color'] = '#B6E880'

        return sample_cols_json, meta_cols_json, col_table.to_dict('records'),\
            button_style

@app.callback(
    Output('renamed_sample_cols', 'data'),
    Output('re_warning', 'children'),
    Output('preview_rename', 'style'),
    Input('preview_rename', 'n_clicks'),
    State('sample_cols', 'children'),
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
        button_style['background-color'] = '#B6E880'

        return col_table.to_dict('records'), '', button_style

@app.callback(
    Output('processed_table', 'children'),
    Output('process_table_button', 'style'),
    Input('process_table_button', 'n_clicks'),
    State('raw_table', 'children'),
    State('sample_cols', 'children'),
    State('meta_cols', 'children'),
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
    State('process_table_button', 'style'),
)

def process_table(n_clicks, raw_table_json, sample_cols_json,
    meta_cols_json, search_re, replacement_re, filter_rows,
    rename_samples, log_transform, transform_opt,
    remove_incomplete_rows, merge_reps, merge_opt, replace_nulls,
    replace_opt, impute_dist, impute_width, button_style):
    
    
    if n_clicks is None:
        raise PreventUpdate
    else:
        # parse raw table and columns from json formats
        raw_table = pd.read_json(raw_table_json)
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
            proteingroup=raw_table
        )
        # designate raw table as filtered_table in the class
        # for specific exception to this Dash App
        # and use it as a default output table
        ms_tables.filtered_table = raw_table
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
        
        # convert final table into json to store in frontend
        final_table_json = final_table.to_json()

        # change button color when finished
        button_style['background-color'] = '#B6E880'

        return final_table_json, button_style

@app.callback(
    Output('download-ms-table-csv', 'data'),
    Output('download_table', 'style'),
    Input('download_table', 'n_clicks'),
    State('processed_table', 'children'),
    State('ms_save_name', 'value'),
    State('download_table', 'style'),
    prevent_initial_call=True
)
def download_ms_table(n_clicks, table_json, save_name, button_style):
    """
    load table from the frontend json format, and save in a csv
    """
    download = pd.read_json(table_json)

    # json reads multi-index tuples as literal strings
    # so convert back to proper multi-level columns
    column_tuples = [eval(name) for name in list(download)]
    download.columns = pd.MultiIndex.from_tuples(column_tuples)


    button_style['background-color'] = '#B6E880'

    return dcc.send_data_frame(download.to_csv,
        save_name), button_style

@app.callback(
    Output('download-configs-csv', 'data'),
    Output('download_configs', 'style'),
    Input('download_configs', 'n_clicks'),
    State('raw_table_sep', 'value'),
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
    State('download_configs', 'style'),
    State('configs_save_name', 'value'),
    prevent_initial_call=True
)
def download_configs(n_clicks, table_sep, top_row_skip, sample_search_re,
    search_re, replacement_re, filter_rows,
    rename_samples, log_transform,
    transform_opt, remove_incomplete_rows, merge_reps, merge_opt,
    replace_nulls, replace_opt, impute_dist, impute_width,
    button_style, save_name):

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

    button_style['background-color'] = '#B6E880'


    return dcc.send_data_frame(configs.to_csv,
        save_name), button_style

@app.callback(
    Output('raw_table_sep', 'value'),
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
    Output('configs_save_name', 'value'),
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
        
        # eval dict for interpreting strings as lists in read_csv
        converters = {
            'filter_rows': pd.eval,
            'rename_samples': pd.eval,
            'log_transform': pd.eval,
            'merge_replicates': pd.eval,
            
        }

        decoded = base64.b64decode(content_string)


        ms_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')))


if __name__ == "__main__":
    app.run_server(debug=True)