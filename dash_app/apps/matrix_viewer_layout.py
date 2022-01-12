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
        html.P('Clustergram viewer',
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
                    id='mat_raw_table_upload',
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
                # html.P('Uploaded feature table', id='mat_raw_table_filename', 
                #     style={
                #         'textAlign': 'center',
                #         'border' : '0.5px #BDC3C7 solid',
                #         'overflowX': 'auto',
                #         'width': '60%',
                #         'marginLeft': '20%',
                #         'marginTop': '1%'
                #     }
                # ),
                
                html.Div([
                    html.Button('Read table!', id = 'mat_read_table_button')],
                    style={
                        'marginTop': '2%',
                        'marginLeft': '30%', 
                        'width': '40%',
                        'verticalAlign': 'top'
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

        
        html.Div([
            html.Div([
                html.P('Select a label', style={'textAlign': 'center',
                    'fontSize':18}),
            ],
                style={
                    'vertical-align': 'top',
                    'marginLeft': '10%',
                    'width': '80%' ,
                }
            ),
            dcc.Dropdown(id='label_select',
                placeholder='Select a label',
                style={
                    'textAlign': 'left',
                    'width': '90%',
                    'marginLeft': '5%'
                }
            ),
            html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),

            html.P('Select columns for the the matrix', style={'textAlign': 'center',
                'fontSize': 18}),
            dcc.Checklist(
                id='features_checklist',
                style={
                    'overflowY': 'auto',
                    'overflowX': 'auto',
                    'height': '250px',
                    'marginLeft': '10%',
                    'width': '80%',
                    'border': '0.5px #BDC3C7 solid',
                }),
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '390px',
                'marginTop': '1%',
                'width': '29%',
                'borderRight': '1px #A6ACAF solid'
            }),
        

        html.Div([
            html.Div([
                html.P('Feature scaling',
                    style={'textAlign': 'center',
                        'fontSize': 18,
                        'marginTop':'0%',
                        'marginBottom': '0%'}
                ),
                dcc.Dropdown(id='feature_scaling',
                    options=[
                        {'label': 'None', 'value':'None'},
                        {'label': 'Standardization', 'value':'standard'},
                        {'label': 'Robust Scaler', 'value':'robust'},
                        {'label': 'Min-max scaling', 'value':'min_max'},
                        {'label': 'l1 normalization', 'value':'l1_norm'},
                        {'label': 'l2 normalization', 'value':'l2_norm'},
                    ],
                    value='standard',
                    style={
                        'marginLeft': '10%',
                        'textAlign': 'center',
                        'width': '80%'
                    }
                )
            ], style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'width': '50%'

            }),
            html.Div([
                html.Button('Scale Data', id = 'scale_data_button',
                    style={
                    'width':'70%',
                    'marginTop':'5%'
                    }
                )


            ], style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'width': '50%'
            }),

            html.Hr(style={'marginTop':'2%', 'marginBottom':'2%'}),

            html.P('Summary Statistics',
                    style={'textAlign': 'center',
                        'fontSize': 18}
                ),
            dash_table.DataTable(
                id='data_metrics',
                columns=[
                    {'name':'Minimum value', 'id': 'min'},
                    {'name':'Maximum value', 'id': 'max'},
                    {'name':'Average', 'id': 'avg'},
                    {'name':'Stdev', 'id': 'stdev'},
                    ],
                style_cell_conditional=[
                    {'if': {'column_id': 'min'},
                    'width': '25%'},
                    {'if': {'column_id': 'max'},
                    'width': '25%'},
                    {'if': {'column_id': 'avg'},
                    'width': '25%'},
                    {'if': {'column_id': 'stdev'},
                    'width': '25%'},
                ],
                data=[{
                    'min': [0],
                    'max': [0],
                    'avg': [0],
                    'stdev': [0]
                }],
                style_header={
                    'fontWeight':'bold'
                },
                style_cell={
                    'textAlign': 'center'
                },
                style_table={
                    'marginLeft': '10%',
                    'height': '80px',
                    'width': '80%'
                },
                fixed_rows={'headers':True, 'data':0}),

            html.Hr(style={'marginTop':'2%', 'marginBottom':'2%'}),

            html.P('Customize Color Scale',
                    style={'textAlign': 'center',
                        'marginBottom': '0.5%',
                        'fontSize': 18}
                ),
            html.Div([
                html.P('Minimum Value',
                    style={'textAlign': 'center',
                        'marginBottom': '0.5%',
                        'fontSize': 15}
                ),
                dcc.Input(
                    id = 'colorscale_min',
                    type = 'number',
                    value=0,
                    style = {'width': '80%',
                        'marginLeft': '10%'}
                ),

            ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'width': '22.5%'
            }),
            html.Div([
                html.P('Maximum Value',
                    style={'textAlign': 'center',
                        'marginBottom': '0.5%',
                        'fontSize': 15}
                ),
                dcc.Input(
                    id = 'colorscale_max',
                    type = 'number',
                    value=100,
                    style = {'width': '80%',
                        'marginLeft': '10%'}
                ),
            ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'width': '22.5%'
            }),
            html.Div([
                html.P('Color map',
                    style={'textAlign': 'center',
                        'marginBottom': '0.5%',
                        'fontSize': 15}
                ),
                dcc.Dropdown(id='colormap',
                    options=[
                        {'label': 'Perseus', 'value':'perseus'},
                        {'label': 'Viridis', 'value':'Viridis'},
                        {'label': 'Cividis', 'value':'Cividis'},
                        {'label': 'Plasma', 'value': 'Plasma'},
                        {'label': 'Blues', 'value': 'Blues'},
                        {'label': 'YellowGreen', 'value': 'YlGn'}
                    ],
                    value='perseus',
                    style={
                        'marginLeft': '5%',
                        'textAlign': 'left',
                        'width': '90%'
                    }
                )

            ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'width': '22.5%'
            }),
            html.Div([
                html.Button('Apply options!', id='color_button',
                style={'width': '80%',
                    'marginLeft':'10%',
                    'marginTop': '7%'})
            ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'width': '32.5%'
            }),
            dcc.Graph(id='color_bar', style={
                'height': '40px',
                'width':'80%',
                'marginLeft':'10%',
                'marginTop':'2%'})


        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '390px',
                'width': '70%'
            }
        ),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.Div([
            html.P('Clustergram Options',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginTop':'15%'}
            ),
            html.Hr(style={'marginTop':'8%', 'marginBottom':'8%'}),

            html.P('Sample clustering',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginBottom': '1%'}
            ),
            html.P('Observations are automatically clustered',
                style={'textAlign': 'center',
                    'fontSize': 13,
                    'marginBottom': '4%'}
            ),

            dcc.Checklist(
                id='cluster_checks',
                options=[
                    {'label': 'Cluster individual samples/replicates', 'value': 'bait_clust'},
                    {'label': 'Cluster samples grouped by replicates', 'value': 'group_clust'},
                ],
                value=[]
            ),
            
            html.Hr(style={'marginTop':'8%', 'marginBottom':'8%'}),

            html.P('Tick labels',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginBottom': '4%'}
            ),

            dcc.Checklist(
                id='tick_checks',
                options=[
                    {'label': 'Show sample labels', 'value': 'sample_ticks'},
                    {'label': 'Show obs. labels', 'value': 'obs_ticks'},
                ],
                value=['sample_ticks', 'obs_ticks']
            ),


            html.Hr(style={'marginTop':'8%', 'marginBottom':'8%'}),
            html.Div([
                html.Button('Create plot!', id = 'generate_matrix',
                    style={
                    'width':'95%',
                    'marginLeft': '2.5%'
                    }
                )
            ],
                style={
                    'marginTop': '2%',
                    'marginLeft': '5%', 
                    'width': '90%',
                    'verticalAlign': 'top'
                }),

        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '700px',
                'width': '15%',
                'borderRight': '1px #A6ACAF solid',
            }
        ),
        html.Div([
            dcc.Graph(id='matrix_fig', style=
            {'height': '100%'})
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '700px',
                'width': '84%',
            }),
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        # Hiddn divs inside the app for computations and storing intermediate values
        html.Div(
            id='mat_processed_table', style={'display': 'none'}),
 



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