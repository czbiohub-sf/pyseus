import base64
import datetime
from dash.html.Label import Label
import markdown
import io
import json
import simplejson
import dash
from dash.dependencies import Input, Output, State, ALL, MATCH
from dash.exceptions import PreventUpdate
from dash import dcc
from dash import html
from dash import dash_table

import pandas as pd
import numpy as np

from dapp import app
from dapp import saved_processed_table, cycle_style_colors



# This is the data set-up for the table of "processed_table_status"
empty_table = pd.DataFrame()
slots = []
tables = []
for i in np.arange(1, 7):
    slots.append(i)
    tables.append('No table processed')
empty_table['slots'] = slots
empty_table['table_name'] = tables
empty_table_dict = empty_table.to_dict('records')


layout = html.Div([
    html.P('the Pyseus Explorer',
        style={'textAlign': 'center', 'fontSize': 28, 'marginTop': '2%',
            'marginBottom': '1%'}),
    html.Hr(),
    html.P('This web-app package includes data exploratory tools specific for mass-spec outputs.\
        Please use the drop-down menu on the top-right to navigate to different tools.\
        For documentation on how to use the tools, please visit here:',
        style={'textAlign': 'center', 'fontSize': 18, 'marginTop': '2%',
            'marginBottom': '1%', 'marginLeft': '15%', 'width': '70%'}),
    dcc.Link('Pyseus Read-the-docs', href='https://pyseus-czb.readthedocs.io/en/main/',
        style={'textAlign': 'center', 'fontSize': 18,
            'marginBottom': '1%', 'marginLeft': '65%', 'width': '35%'}),

    html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),


    html.Div([

        html.P('Processed Tables Status',
            style={'textAlign': 'center', 'fontSize': 20,
                'marginBottom': '1%', 'marginLeft': '20%', 'width': '60%'}),
        dash_table.DataTable(
            id='processed_table_status',
            columns=[
                {'name': 'Slots', 'id': 'slots'},
                {'name': 'Table name', 'id': 'table_name'}],

            style_cell_conditional=[
                {'if': {'column_id': 'slots'},
                'width': '20%'},
            ],
            data=empty_table_dict,
            style_header={
                'fontWeight': 'bold'
            },
            style_cell={
                'textAlign': 'center'
            },
            style_table={
                'marginLeft': '10%',
                'overflowX': 'auto',
                'height': '210px',
                'width': '80%',
                'border' : '0.5px #BDC3C7 solid'
            },
            fixed_rows={'headers': True, 'data': 0})],
        style={
            'display': 'inline-block', 'width': '49%', 'borderRight': '1px #A6ACAF solid',
            'verticalAlign': 'Top', 'height': '260px'}

    ),

    html.Div([
        html.P('Upload a processed table',
            style={'textAlign': 'center', 'fontSize': 20, 'marginTop': '1%',
                'marginBottom': '1%', 'marginLeft': '20%', 'width': '60%'}),

        dcc.RadioItems(id='download_slot', options=[
            {'label': 'Slot 1', 'value': 0},
            {'label': 'Slot 2', 'value': 1},
            {'label': 'Slot 3', 'value': 2},
        ], labelStyle={'display': 'inline-block'},
            value=0,
            style={'marginLeft': '35%', 'width': '60%', 'marginTop': '2%', 'marginBottom': '1%',
                'fontSize': 16}),

        dcc.Download(id='home_dl'),
        dcc.Upload(
            id='home_raw_table_upload',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select a standard table')
            ]),
            style={
                'marginLeft': '15%',
                'width': '70%',
                'height': '70px',
                'lineHeight': '70px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'background-color': '#E5E7E9'
            },
            # Only single file
            multiple=False
        ),
        html.Button('Upload Table', id='home_ul_button', style={
            'marginLeft': '15%', 'width': '70%', 'marginTop': '4%',
            'height': '5%'
        }),




    ], style={
        'display': 'inline-block', 'width': '50%', 'verticalAlign': 'Top'}
    ),


    html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

    html.P('Developed by Kibeom Kim, for internal use at CZ Biohub only.',
        style={'textAlign': 'center', 'fontSize': 16, 'marginTop': '2%',
            'marginBottom': '1%', 'marginLeft': '15%', 'width': '70%'}),


])


@app.callback(
    Output('home_raw_table_upload', 'children'),
    Output('home_raw_table_upload', 'style'),
    Input('home_raw_table_upload', 'filename'),
    State('home_raw_table_upload', 'style')
)
def display_upload_ms_filename(filename, style):
    if filename is None:
        raise PreventUpdate
    else:
        style = cycle_style_colors(style)
        return filename, style


@app.callback(

    Output('slot_label_1', 'children'),
    Output('slot_label_2', 'children'),
    Output('slot_label_3', 'children'),
    Output('home_ul_button', 'style'),

    Input('home_ul_button', 'n_clicks'),

    # States from home page upload
    State('home_raw_table_upload', 'contents'),
    State('home_raw_table_upload', 'filename'),
    State('download_slot', 'value'),
    State('session_id', 'data'),


    State('slot_label_1', 'children'),
    State('slot_label_2', 'children'),
    State('slot_label_3', 'children'),
    State('home_ul_button', 'style'),

    prevent_initial_call=True
)
def upload_table(n_clicks, content,
        filename, slot_num, session_id, label_1, label_2, label_3,
        home_button_style):
    """
    When a user uploads a table and clicks upload button, save the table
    and table name to the server under the right slot
    """

    if n_clicks is None:
        raise PreventUpdate

    # unique session id for the specific slot the table will be cached
    session_slot = session_id + str(slot_num)

    # assign labels to a list for easy slot designation
    labels = [label_1, label_2, label_3]


    # parse file
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)

    raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
        low_memory=False, header=[0, 1], index_col=0)

    # cache the table to the specific slot
    _ = saved_processed_table(session_slot, raw_table, overwrite=True)

    # Save the filename to the specific slot
    labels[slot_num] = filename

    home_button_style = cycle_style_colors(home_button_style)

    return labels[0], labels[1], labels[2], home_button_style


@app.callback(
    Output('processed_table_status', 'data'),
    Input('slot_label_1', 'children'),
    Input('slot_label_2', 'children'),
    Input('slot_label_3', 'children'),
    Input('slot_label_4', 'children'),
    Input('slot_label_5', 'children'),
    Input('slot_label_6', 'children'),
    State('processed_table_status', 'data'),
)
def update_table_status(label_1, label_2, label_3, label_4, label_5, label_6, data):
    """
    Automatic update of file/table names in the upload/save slots.
    """
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
