import base64
import datetime
import markdown
import io
import json
import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dcc
from dash import html
from dash import dash_table

import pandas as pd
import plotly.express as px
import sys

from preprocessing_layout import create_layout

sys.path.append('../../')
from pyseus import basic_processing as bp


# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# App Layout
app.layout = create_layout()

# @app.callback(
#     Output('upload-status', 'children'),
#     Input('pGroup-upload', 'contents'),
#     State('pGroup-upload', 'filename')
#     )
# def check_filename(contents, filename):
#     """
#     Check if the file has 'proteinGroups' or '.txt' in the filename
#     """
#     if contents is not None:
#         if 'proteinGroups' not in filename:
#             return "File name does not include 'proteinGroups', check the file again"
#         elif 'txt' not in filename:
#             return "File is not in '.txt' format, check the file again"
#         else:
#             return "Upload file OK!"

@app.callback(
    Output('raw_table', 'children'),
    Output('all_cols', 'children'),
    Output('sample_cols_checklist', 'options'),
    Output('meta_cols_checklist', 'options'),
    Input('read_table_button', 'n_clicks'),
    State('raw_table_upload', 'contents'),
    State('raw_table_sep', 'value'),
    State('skip_top_rows', 'value')
    )
def parse_raw_table(n_clicks, content, sep_type, num_skip_row):
    """
    initiate QualityControl class with the uploaded proteingroups file
    """

    if n_clicks is None:
        raise PreventUpdate
    else:
        # parse txt (tsv) file as pd df from upload
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        pgroups = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t',
            low_memory=False)
        
        # initiate RawTable class
        exp = bp.RawTables(
            experiment_dir='',
            pg_file='',
            intensity_type='Intensity ',
            proteingroup=pgroups,
            file_designated=True)
        exp.filter_table()
        filtered_table = exp.filtered_table.to_json()
        intensity_cols = json.dumps(exp.intensity_cols)
        info_cols = json.dumps(exp.info_cols)

        return filtered_table, intensity_cols, info_cols

# @app.callback(
#     Output('original_intensity_cols', 'data'),
#     Input('intensity_cols', 'children'),
#     Input('renamed_cols', 'children'),
#     State('button2', 'n_clicks'),
# )
# def convert_intensity_cols_to_markdown(orig_cols_json, renamed_cols_json,
#     n_clicks):
#     """
#     from the loaded QualityControl class, extract 
#     the list of intensity cols and convert to markdown
#     """
#     if orig_cols_json is None:
#         raise PreventUpdate
#     else:
#         if renamed_cols_json is None:
#             intensity_cols = json.loads(orig_cols_json)
#         else:
#             intensity_cols = json.loads(renamed_cols_json)
#         col_table = pd.DataFrame()
#         col_table['column names'] = intensity_cols

    
#         return col_table.to_dict('records')


# @app.callback(
#     Output('renamed_intensity_cols', 'data'),
#     Output('re_warning', 'children'),
#     Input('button1', 'n_clicks'),
#     State('intensity_cols', 'children'),
#     State('search_re', 'value'),
#     State('replacement_re', 'value')
# )
# def test_rename_cols(n_clicks, cols_json, search_re, replacement_re):
#     if n_clicks is None:
#         raise PreventUpdate

#     else:
#         search_re = search_re.split(';')
#         replacement_re = replacement_re.split(';')
#         if len(search_re) != len(replacement_re):
#             return None, "Number of search and replacement RE(s) does not match!"
#         replacement_re = ['' if x == 'NONE' else x for x in replacement_re ]
#         intensity_cols = json.loads(cols_json)
#         new_cols = bp.sample_rename(intensity_cols,
#             search_re, replacement_re)
        
#         col_table = pd.DataFrame()
#         col_table['new column names'] = new_cols

#         return col_table.to_dict('records'), ''


# @app.callback(
#     Output('renamed_table', 'children'),
#     Output('renamed_cols', 'children'),
#     Input('button2', 'n_clicks'),
#     State('filtered_table', 'children'),
#     State('intensity_cols', 'children'),
#     State('search_re', 'value'),
#     State('replacement_re', 'value')
#     )
# def rename_cols(n_clicks, filtered_table_json, cols_json,
#     search_re, replacement_re):
#     if n_clicks is None:
#         raise PreventUpdate
#     else:
#         search_re = search_re.split(';')
#         replacement_re = replacement_re.split(';')
#         replacement_re = ['' if x == 'NONE' else x for x in replacement_re ]
        
#         intensity_cols = json.loads(cols_json)
#         filtered_table = pd.read_json(filtered_table_json)

#         # initiate RawTables class
#         exp = bp.RawTables(
#             experiment_dir='',
#             pg_file='',
#             intensity_type='Intensity ',
#             proteingroup=None,
#             file_designated=True)
        
#         exp.filtered_table = filtered_table
#         exp.intensity_cols = intensity_cols
#         exp.rename_columns(search_re, replacement_re)

#         renamed_table = exp.filtered_table.to_json()
#         new_cols = json.dumps(exp.intensity_cols)

#         return renamed_table, new_cols





if __name__ == "__main__":
    app.run_server(debug=True)