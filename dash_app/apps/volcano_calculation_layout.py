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
from numpy import true_divide

import pandas as pd
import plotly.express as px

import sys
sys.path.append('../../')

def calculation_layout():
    return html.Div([

        html.Div([
            
            html.P('Upload a processed feature table or use a pre-loaded one',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom':'0%',
                    'marginTop': '4%'}),
            html.P("""The script requires a standard format
                of feature table. Please use the Process Table page to create one.""",
                style={
                    'textAlign': 'center',
                    'fontSize': 14,
                    'marginTop': '0%'}),

            # upload component
            html.Div([
                dcc.Upload(
                    id='vol_raw_table_upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select a feature table')
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
                    html.Button('Read table!', id = 'vol_read_table_button')],
                    style={
                        'marginTop': '2%',
                        'marginLeft': '30%', 
                        'width': '40%',
                        'verticalAlign': 'top',
                        'white-space': 'normal'
                    }),
            ],
            style = {
                'display': 'inline-block',
                'borderRight': '1px #A6ACAF solid',
                'width': '49%',
                'textAlign': 'center',
                'vertical-align': 'top',
                }
            ),
            html.Div([

                dcc.Dropdown(id='vol_preloaded_dropdown',
                    options=[
                        {'label': 'slot 1', 'value': '1'},
                        {'label': 'slot 2', 'value': '2'},
                        {'label': 'slot 3', 'value': '3'},
                    ],
                    placeholder='Select a pre-loaded table',
                    style={
                        'textAlign': 'center',
                        'width': '90%',
                        'marginLeft': '5%',
                        'marginTop': '2%'
                    }
                ),
                html.Div([
                    html.Button(
                        'Load data!',
                        id =  'vol_preload_button',)
                    ], style={
                        'marginTop': '2%',
                        'textAlign': 'center',
                        'display': 'inline-block',
                        'width': '40%',
                        'marginLeft': '30%'

                    }),


            ],
            style = {
                'display': 'inline-block',
                'width': '50%',
                'vertical-align': 'top',
                }
            ),
        ]),
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.P('Enrichment / Significance calculation',
                style={'textAlign': 'center',
                    'fontSize':20, 'marginBottom':'1.5%'}),
        html.Div([
            html.Div([
                html.P('Select a control selection mode',
                    style={'textAlign': 'center',
                        'fontSize':18,
                        'marginBottom':'0%'}),
                html.P('By default, manual control selection includes all samples \
                    except the sample itself.',
                    style={'textAlign': 'center',
                        'fontSize':14,
                        'marginTop':'0%',
                        'lineHeight':'15px'})
            ],
                style={
                    'vertical-align': 'top',
                    'marginLeft': '10%',
                    'width': '80%' ,
                }
            ),

            dcc.RadioItems(id='calculation_options',
                options=[
                    {'label': 'Manually select control samples',
                    'value': 'manual'},
                    {'label': 'Automatically calulate the null model',
                    'value': 'automatic'}
                    ],
        
                style={
                    'textAlign': 'left',
                    'width': '90%',
                    'marginLeft': '5%'
                },
                value='manual'
            ),


            html.Hr(style={'marginTop': '4%', 'marginBottom': '4%'}),
            html.P('Download current control matrix', style={'textAlign': 'center',
                'fontSize': 18}),
            html.Button('Download', id='download_matrix_button',
                style={
                    'marginLeft':'7.5%',
                    'width': '85%',
                    'white-space': 'normal'
                }
            ),
            dcc.Download(id="download_control_matrix"),

            html.Hr(style={'marginTop': '4%', 'marginBottom': '4%'}),
            html.P('Import a custom control matrix', style={'textAlign': 'center',
                'fontSize': 18}),
            dcc.Upload(
                id='null_matrix_upload',
                children=html.Div([
                    'D&D or ',
                    html.A('Select a control sample matrix')
                ]),
                style={
                    'marginLeft': '7.5%',
                    'width': '85%',
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
            html.Button('Apply uploaded matrix', id='null_matrix_button',
                style={
                    'marginLeft':'7.5%',
                    'width': '85%',
                    'marginTop': '4%',
                    'white-space': 'normal'
                }
            ),

        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '420px',
                'marginTop': '1%',
                'width': '32%',
                'borderRight': '1px #A6ACAF solid'
            }),
        

        html.Div([
            html.P('Review controls for a selected sample', style={'textAlign': 'center',
                'fontSize': 18, 'marginTop': '4%'}),
            
            dcc.Dropdown(id='view_sample_controls',
                placeholder='Select a sample',
                style={
                    'textAlign': 'center',
                    'width': '90%',
                    'marginLeft': '5%',
                    'marginTop': '2%',
                    'marginBottom':'3%'
                }
            ),
            dash_table.DataTable(
                id='sample_controls',
                columns=[{'name':'control sample names', 'id': 'sample_control'}],
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
                        'height': '300px',
                        'width': '80%',
                        'border' : '0.5px #BDC3C7 solid'
                    },
                    fixed_rows={'headers':True, 'data':0}
                ),
            ],

            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '420px',
                'width': '32%',
                'borderRight': '1px #A6ACAF solid'
            }
        ),
        html.Div([
            html.P('Samples to edit', style={'textAlign': 'center',
                'fontSize': 16, 'lineHeight':'10px'}),
            dcc.Checklist(
                id='edit_samples',
                style={
                    'overflowY': 'auto',
                    'overflowX': 'auto',
                    'height': '140px',
                    'marginLeft': '10%',
                    'width': '80%',
                    'border': '0.5px #BDC3C7 solid',
                }),
            html.P('Select controls for the chosen samples', style={'textAlign': 'center',
                'fontSize': 16, 'lineHeight':'15px', 'marginTop':'2%'}),
            dcc.Checklist(
                id='edit_controls',
                style={
                    'overflowY': 'auto',
                    'overflowX': 'auto',
                    'height': '170px',
                    'marginLeft': '10%',
                    'marginTop':'1%',
                    'width': '80%',
                    'border': '0.5px #BDC3C7 solid',
                }),
            html.Button('Apply control selection', id='control_apply_button',
                style={
                    'marginLeft':'7.5%',
                    'width': '85%',
                    'marginTop': '4%',
                    'white-space': 'normal'
                }
            )
           
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '420px',
                'width': '33%'
            }
        ),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.Div([

            html.P('Enrichment option', style={'textAlign': 'center',
                'fontSize': 18, 'lineHeight':'15px', 'marginTop':'2%'}),

            dcc.RadioItems(id='enrichment_option',
                options=[
                    {'label': 'Absolute enrichment',
                    'value': 'False'},
                    {'label': 'Relative enrichment (stdev units)',
                    'value': 'True'}
                    ],
        
                style={
                    'textAlign': 'center',
                    'width': '90%',
                    'marginLeft': '5%'
                },
                value='True'
            ),
            ],

            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '100px',
                'width': '50%',
            }
        ),

        html.Div([

            html.Button('Calculate enrichment!', id='calculate_button',
                style={
                    'marginLeft':'0%',
                    'width': '85%',
                    'marginTop': '0%',
                    'white-space': 'normal',
                    'background-color': 'white'
                }
            ),
            dcc.Download(id="download_pval_table"),
            html.Button('Download calculated table!', id='download_pval_button',
                style={
                    'marginLeft':'0%',
                    'width': '85%',
                    'marginTop': '2%',
                    'white-space': 'normal'
                }
            )

            ],

            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '100px',
                'width': '50%',
            }
        ),





        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        # Hiddn divs inside the app for computations and storing intermediate values
        html.Div(
            id='vol_processed_table', style={'display': 'none'}),
        html.Div(
            id='control_matrix', style={'display': 'none'}
        ),
        html.Div(
            id='samples', style={'display': 'none'}
        ),
        html.Div(
            id='significance_table', style={'display': 'none'})



        # html.P('Uploaded feature table', id='raw_table_filename', 
        #             style={
        #                 'textAlign': 'center',
        #                 'border' : '0.5px #BDC3C7 solid',
        #                 'overflowX': 'auto',
        #                 'width': '60%',
        #                 'marginLeft': '20%',
        #                 'marginTop': '1%'
        #             }
        #         ),


    ])