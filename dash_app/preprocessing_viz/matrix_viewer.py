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

from matrix_viewer_layout import create_layout

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

# @app.callback(
#     Output('raw_table_upload', 'children'),
#     Output('raw_table_upload', 'style'),
#     Input('raw_table_upload', 'filename'),
#     State('raw_table_upload', 'style')
# )
# def display_upload_ms_filename(filename, style):
#     if filename is None:
#         raise PreventUpdate
#     else:
#         style['background-color'] = '#B6E880'
#         return filename, style


# @app.callback(
#     Output('processed_table', 'children'),
#     Output('features', 'children'),
#     Output('annots', 'children'),
#     Output('table_dims', 'children'),
#     Output('features_checklist', 'options'),
#     Output('features_checklist', 'value'),
#     Output('label_select', 'options'),
#     Output('label_select', 'value'),
#     Output('read_table_button', 'style'),
#     Input('read_table_button', 'n_clicks'),
#     State('raw_table_upload', 'contents'),
#     State('read_table_button', 'style'),
#     )
# def parse_raw_table(n_clicks, content, button_style):
#     """
#     initiate QualityControl class with the uploaded proteingroups file
#     """

#     if n_clicks is None:
#         raise PreventUpdate
#     else:
#         # parse txt (tsv) file as pd df from upload
#         content_type, content_string = content.split(',')
#         decoded = base64.b64decode(content_string)

#         raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
#             low_memory=False, header=[0,1], index_col=0)
        

#         # regular table, features, and annots for UMAP
#         processed_table = raw_table.droplevel(level=0, axis=1)
#         features = list(raw_table['sample'])
#         labels = list(raw_table['metadata'])
#         dims = list(raw_table['sample'].shape)

#         features_json = json.dumps(features)
#         labels_json = json.dumps(labels)
#         dims_json = json.dumps(dims)


#         # feature checklist options 
#         features_opts = [{'label': feature, 'value': feature}
#             for feature in features]


#         # labels/annots options
#         label_opts = []
#         label_opts.append({'label': 'None', 'value': 'None'})
#         for label in labels:
#             label_opts.append({'label': label, 'value': label})

#         annot_val = 'None'

#         # drop duplicates
#         processed_table.drop_duplicates(inplace=True)
#         processed_table_json = processed_table.to_json()


#         if button_style is None:
#             button_style = {}
#         button_style['background-color'] = '#B6E880'

#         return processed_table_json, features_json, labels_json,\
#             dims_json, button_style




    


if __name__ == "__main__":
    app.run_server(debug=True)
