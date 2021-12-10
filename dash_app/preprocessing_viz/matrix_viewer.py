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
from pyseus.plotting import plotly_heatmap as ph

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
    Output('features_checklist', 'options'),
    Output('features_checklist', 'value'),
    Output('label_select', 'options'),
    Output('data_metrics', 'data'),
    Output('color_button', 'n_clicks'),
    Output('read_table_button', 'style'),
    Input('read_table_button', 'n_clicks'),
    State('raw_table_upload', 'contents'),
    State('read_table_button', 'style'),
    State('color_button', 'n_clicks'),
    prevent_initial_call=True


    )
def parse_raw_table(n_clicks, content, button_style, color_clicks):
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
        labels = list(raw_table['metadata'])


        # feature checklist options 
        features_opts = [{'label': feature, 'value': feature}
            for feature in features]


        # labels/annots options
        label_opts = []
        label_opts.append({'label': 'None', 'value': 'None'})
        for label in labels:
            label_opts.append({'label': label, 'value': label})

        processed_table_json = processed_table.to_json()

        if button_style is None:
            button_style = {}
        button_style['background-color'] = '#B6E880'

        if color_clicks is None:
            color_clicks = 1
        else:
            color_clicks += 1

        matrix = processed_table[features]
        metrics = [{
            'min': [np.round(matrix.values.min(),2)],
            'max': [np.round(matrix.values.max(),2)],
            'avg': [np.round(matrix.values.mean(),2)],
            'stdev': [np.round(matrix.values.std(),2)]
        }]

        return processed_table_json, features_opts, features,\
            label_opts, metrics, color_clicks, button_style



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
    


if __name__ == "__main__":
    app.run_server(debug=True)
