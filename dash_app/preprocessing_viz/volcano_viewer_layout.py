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
        html.P('Volcano plot generator',
            style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'2%',
                'marginBottom': '1%'}),
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.Div([
            
            html.P('Upload a processed feature table or use a pre-loaded one',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom':'0%'}),
            html.P("""The UMAP generating script requires a standard format
                of feature table. Please use the Process Table page to create one.""",
                style={
                    'textAlign': 'center',
                    'fontSize': 14,
                    'marginTop': '0%'}),

            # upload component
            html.Div([
                dcc.Upload(
                    id='raw_table_upload',
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
                    html.Button('Read table!', id = 'read_table_button')],
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

                dcc.Dropdown(id='preloaded_dropdown',
                    options=[
                        {'label': 'table 1', 'value': 'table1'},
                        {'label': 'table 2', 'value': 'table2'},
                        {'label': 'table 3', 'value': 'table3'},
                        {'label': 'table 4', 'value': 'table4'},
                        {'label': 'table 5', 'value': 'table5'},
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
                        id =  'load-preload-button',)
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


            html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),

            html.P('Import a manual control sample matrix', style={'textAlign': 'center',
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
            html.Button('Apply matrix', id='null_matrix_button',
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
                'marginTop': '1%',
                'width': '32%',
                'borderRight': '1px #A6ACAF solid'
            }),
        

        html.Div([
            html.P('View selected controls', style={'textAlign': 'center',
                'fontSize': 18}),
            
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
            html.P('Edit controls for selected samples', style={'textAlign': 'center',
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



        # html.Div([
        #     html.P('Clustergram Options',
        #         style={'textAlign': 'center',
        #             'fontSize': 20,
        #             'marginTop':'15%'}
        #     ),
        #     html.Hr(style={'marginTop':'8%', 'marginBottom':'8%'}),

        #     html.P('Sample clustering',
        #         style={'textAlign': 'center',
        #             'fontSize': 15,
        #             'marginBottom': '1%'}
        #     ),
        #     html.P('Observations are automatically clustered',
        #         style={'textAlign': 'center',
        #             'fontSize': 13,
        #             'marginBottom': '4%'}
        #     ),

        #     dcc.Checklist(
        #         id='cluster_checks',
        #         options=[
        #             {'label': 'Cluster individual samples/replicates', 'value': 'bait_clust'},
        #             {'label': 'Cluster samples grouped by replicates', 'value': 'group_clust'},
        #         ],
        #         value=[]
        #     ),
            
        #     html.Hr(style={'marginTop':'8%', 'marginBottom':'8%'}),

        #     html.P('Tick labels',
        #         style={'textAlign': 'center',
        #             'fontSize': 15,
        #             'marginBottom': '4%'}
        #     ),

        #     dcc.Checklist(
        #         id='tick_checks',
        #         options=[
        #             {'label': 'Show sample labels', 'value': 'sample_ticks'},
        #             {'label': 'Show obs. labels', 'value': 'obs_ticks'},
        #         ],
        #         value=['sample_ticks', 'obs_ticks']
        #     ),


        #     html.Hr(style={'marginTop':'8%', 'marginBottom':'8%'}),
        #     html.Div([
        #         html.Button('Create plot!', id = 'generate_matrix',
        #             style={
        #             'width':'95%',
        #             'marginLeft': '2.5%'
        #             }
        #         )
        #     ],
        #         style={
        #             'marginTop': '2%',
        #             'marginLeft': '5%', 
        #             'width': '90%',
        #             'verticalAlign': 'top'
        #         }),

        # ],
        #     style={
        #         'vertical-align': 'top',
        #         'display': 'inline-block',
        #         'height': '700px',
        #         'width': '15%',
        #         'borderRight': '1px #A6ACAF solid',
        #     }
        # ),
        # html.Div([
        #     dcc.Graph(id='matrix_fig', style=
        #     {'height': '100%'})
        # ],
        #     style={
        #         'vertical-align': 'top',
        #         'display': 'inline-block',
        #         'height': '700px',
        #         'width': '84%',
        #     }),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        # Hiddn divs inside the app for computations and storing intermediate values
        html.Div(
            id='processed_table', style={'display': 'none'}),
 



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