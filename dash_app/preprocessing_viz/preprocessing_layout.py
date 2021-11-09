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
sys.path.append('../../')

def create_layout():
    return html.Div([
        
        # Header tags
        html.P('MS Table Processing',
            style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'3%',
                'marginBottom': '1%'}),


        # Upload config and proteinGroup section
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        
        html.Div([
            
            html.P('Upload MS table and config file (optional)',
                style={'textAlign': 'center', 'fontSize': 22}),

            # upload component
            html.Div([
                dcc.Upload(
                    id='raw_table_upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select an MS table')
                    ]),
                    style={
                        'marginLeft': '10%',
                        'width': '80%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'background-color': '#E5E7E9'
                    },
                    # Only single file
                    multiple=False
                ),
                html.Div([
                    dcc.RadioItems(
                        id='raw_table_sep',
                        options=[
                            {'label': 'comma-separated', 
                            'value': ','},
                            {'label': 'tab-separated',
                            'value': '\t'}
                        ],
                        value=',',
                        ),
                    ],
                    style={
                        'display': 'inline-block',
                        'marginLeft': '2%',
                        'marginTop': '1.5%',
                        'width': '30%'
                    }),
                html.Div([
                    html.P(
                        '# top rows to skip', 
                        style={
                            'fontSize': 14,
                            'textAlign' : 'center',
                            'marginTop' : '1%',
                            'marginBottom': '0.5%'
                        }
                    ),
                    dcc.Input(
                        id = 'skip_top_rows',
                        value=  0,
                        type = 'number', 
                        style = {
                            'height': '3%',
                            'width': '50%'
                        }
                        ),
                    ],
                    style={
                        'display': 'inline-block',
                        'marginTop': '1.5%',
                        'width': '20%'
                    }),
                html.Div([
                    html.Button('Read table!', id = 'read_table_button')],
                    style={
                        'display': 'inline-block',
                        'marginTop': '4%',
                        'width': '30%',
                        'verticalAlign': 'top'
                    }),
            ],
            style = {
                'display': 'inline-block',
                'borderRight': '1px grey solid',
                'width': '49%',
                'textAlign': 'center',
                'vertical-align': 'top',
                }
            ),
            html.Div([
                dcc.Upload(
                    id='config-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select a config file')
                    ]),
                    style={
                        'marginLeft': '10%',
                        'width': '80%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'background-color': '#E5E7E9'
                    },
                    # Only single file
                    multiple=False),
                html.Div([
                    html.Button(
                        'Pre-load configs',
                        id =  'load-config-button',)
                    ], style={
                        'marginTop': '3%',
                        'textAlign': 'center',
                        'display': 'inline-block',
                        'width': '40%',
                        'marginLeft': '10%'

                    }),
                html.Div([
                    html.Button(
                        'Apply configs!',
                        id =  'apply-config-button',
                        ),
                    ], style={
                        'marginTop': '3%',
                        'textAlign': 'center',
                        'display': 'inline-block',
                        'width': '40%',
                        'marginRight': '10%'

                    }),

            ],
            style = {
                'display': 'inline-block',
                'width': '50%',
                'vertical-align': 'top',
                }
            ),
        ]),


        # Column selection
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),    
        html.P('Select sample and metadata columns',
                style={'textAlign': 'center', 'fontSize': 22}),
        
        html.Div([
            html.Div([
                html.P('Sample columns', style={'textAlign': 'center'}),
                dcc.Checklist(
                    id='sample_cols_checklist',
                    style={
                        'marginLeft': '10%',
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '250px',
                        'width': '80%'
                    }),
                html.P('Regular expression search', style = {'textAlign': 'center'}),
                html.Div([
                    dcc.Input(
                        id = 'sample_search_re',
                        type = 'text',
                        placeholder='python regular expression',
                        style = {'width': '90%'}
                    )
                ],
                style= {
                    'display': 'inline-block',
                    'marginLeft': '10%',
                    'width': '50%'
                }),
                html.Div([
                    html.Button('Apply!',
                        id = 'apply_sample_re_button',
                        style={'width': '90%'},
                    )
                ],
                style= {
                    'display': 'inline-block',
                    'width': '30%'
                }),

            ],
            style = {
                'display': 'inline-block',
                'width': '49%',
                'borderRight': '1px grey solid'
                }
            ),
            html.Div([
                html.P('Metadata columns', style={'textAlign': 'center'}),
                dcc.Checklist(
                    id='meta_cols_checklist',
                    style={
                        'marginLeft': '10%',
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '250px',
                        'width': '80%'
                    })],
                style = {
                    'display': 'inline-block',
                    'width': '50%',
                    'vertical-align': 'top'}
            ),    
        ]), 


        html.Div([
            html.Button(
                'Save feature / metadata columns',
                id=  'save-cols-button',
                style = {
                    'marginTop': '3%',
                    'width': '40%'
                }),
            ], style={
                'textAlign': 'center'
            }
        ),



        # Preview Renaming
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.P('Renaming samples preview',
            style={'textAlign': 'center', 'fontSize': 22}),
        html.P("""Sample names should follow the format of 'experiment_sample_rep#',
            for example: CZBSU01_LAMP1_1 or ATL3_2""",
            style={'textAlign':'center'}),
        html.Div([
            html.Div([
                dash_table.DataTable(
                    id='original_intensity_cols',
                    columns=[{'name':'column names', 'id': 'column names'}],
                    style_header={
                        'fontWeight':'bold'
                    },
                    style_cell={
                        'textAlign': 'center'
                    },
                    style_table={
                        'marginLeft': '10%',
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '250px',
                        'width': '80%'
                    },
                    fixed_rows={'headers':True, 'data':0})],
                style = {
                    'display': 'inline-block',
                    'width': '50%',
                    'height': '250px'}
            ),
            html.Div([
                dash_table.DataTable(
                    id='renamed_intensity_cols',
                    columns=[{'name':'new column names', 'id': 'new column names'}],
                    style_header={
                        'fontWeight':'bold'
                    },
                    style_cell={
                        'textAlign': 'center'
                    },
                    style_table={
                        'marginLeft': '10%',
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '250px',
                        'width': '80%'
                    },
                    fixed_rows={'headers':True, 'data':0})],
                style = {
                    'display': 'inline-block',
                    'width': '50%',
                    'height': '250px'}
            ),    ],
            style={'marginBottom':'1.5%'}),
        html.Div([
            html.Div([
                html.P("Search regular expressions, split by semicolons:")   
            ], style={
                'display': 'inline-block',
                'marginLeft': '5%',
                'width': '25%'}),
            html.Div([
                dcc.Input(
                    id = 'search_re',
                    value = 'Intensity ;',
                    type = 'text',
                    style = {'width': '95%'}
                )
            ], style={
                'display': 'inline-block',
                'marginLeft': '2%',
                'width': '60%'}),

        ]),
        html.Div([
            html.Div([
                html.P("Replacement regular expressions, split by semicolons (type NONE for deletion):")   
            ], style={
                'display': 'inline-block',
                'marginLeft': '5%',
                'width': '25%'}),
            html.Div([
                dcc.Input(
                    id = 'replacement_re',
                    value = 'NONE;',
                    type = 'text',
                    style = {'width': '95%'}
                )
            ], style={
                'display': 'inline-block',
                'marginLeft': '2%',
                'width': '60%'}),

        ]),
        html.Div([
            html.Button('Test renaming samples', id = 'button1'),
            html.Div([
            html.Span(
                id='re_warning',
                style={'color':'red'})
            ])
            ],

            style = {
                'textAlign': 'center',
                'marginTop': '2%',
                'display': 'inline-block',
                'width': '50%'
                }),

        html.Div([
            html.Button('Rename samples!', id = 'button2'),
            html.Div([
                html.Span(id='blank')
            ]),
            ],

            style = {
                'textAlign': 'center',
                'marginTop': '2%',
                'display': 'inline-block',
                'width': '50%'
                }),
            

        # Hiddn divs inside the app for computations and storing intermediate values
        html.Div(
            id='raw_table', style={'display': 'none'}),
        html.Div(
            id='all_cols', style={'display': 'none'}),
        html.Div(
            id='sample_cols', style={'display': 'none'}),
        html.Div(
            id='info_cols', style={'display': 'none'}),
        # html.Div(
        #     id='renamed_sample_cols', children=None, style={'display': 'none'}),
        # html.Div(
        #     id='renamed_table', children=None, style={'display': 'none'}),
        # html.Div(
        #     id='grouped_table', children=None, style={'display': 'none'}),
        
    ])
