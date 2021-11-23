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
        html.P('The UMAP generator',
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
                html.P('Uploaded feature table', id='raw_table_filename', 
                    style={
                        'textAlign': 'center',
                        'border' : '0.5px #BDC3C7 solid',
                        'overflowX': 'auto',
                        'width': '60%',
                        'marginLeft': '20%',
                        'marginTop': '1%'
                    }
                ),
                
                html.Div([
                    html.Button('Read table!', id = 'read_table_button')],
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
                        'marginTop': '4%'
                    }
                ),
                html.Div([
                    html.Button(
                        'Load data!',
                        id =  'load-preload-button',)
                    ], style={
                        'marginTop': '4%',
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
        html.P('Feature table dimensions',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom':'0%',
                    'marginTop': '1%'}),
        html.P('0 features X 0 observations, 0 annotations',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginBottom':'0%',
                    'marginTop': '0%'}),
        html.Div([
            html.Button(
                'Transpose matrix',
                id =  'transpose-button',)
            ], style={
                'marginTop': '1%',
                'textAlign': 'center',
                'width': '40%',
                'marginLeft': '30%'

        }),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.P('Select features, label, and annotations',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom':'0%'}),
        
        html.Div([
                html.P('Select features for UMAP generation', style={'textAlign': 'center',
                    'fontSize': 17}),
                dcc.Checklist(
                    id='features_checklist',
                    style={
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '280px',
                        'width': '90%',
                        'border': '0.5px #BDC3C7 solid',
                    }),
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '340px',
                'marginLeft': '5%',
                'marginTop': '1%',
                'width': '34%',
                'borderRight': '1px #A6ACAF solid'
            }),

        html.Div([
            html.Div([
                html.P('Select a label for markers', style={'textAlign': 'center',
                    'fontSize':17}),
            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'marginLeft': '10%',
                    'width': '30%' ,
                    'marginTop': '1.75%'
                }
            ),
            html.Div([
                dcc.Dropdown(id='label_select',
                    placeholder='Select a label',
                    style={
                        'textAlign': 'left',
                        'width': '80%'
                    }
            ),
            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'width': '60%',
                    'marginTop': '1.25%' 
                }
            ),

            html.Hr(style={'marginTop':'2%', 'marginBottom':'0%'}),
            html.P('Select annotations for color mapping',
                style={'textAlign': 'center',
                    'fontSize':17,
                    'marginTop': '2%'}),
            html.Div([
                dcc.Dropdown(id='annot_select',
                    placeholder='Do not map annotations',
                    style={
                        'textAlign': 'center'
                    }
                ),
                ], style={
                    'vertical-align': 'top',
                    'marginLeft': '20%',
                    'width': '60%' 
                }),

            html.P('OR link external annotations', style={
                    'marginTop':'1%',
                    'fontSize':17,
                    'textAlign': 'center'}),
            dcc.Upload(
                    id='annot_table_upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select an annotation table')
                    ]),
                    style={
                        'marginLeft': '10%',
                        'marginTop': '1.5%',
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
            html.P('csv or tsv files accepted. Table header MUST be in top-row.',
                style={
                    'fontSize':13,
                    'textAlign': 'center'
                }),
            html.Div([
                html.P('Merge-key: feature table',
                    style={
                        'textAlign':'center',
                        'fontSize': 14
                    }),
                dcc.Dropdown(id='merge-key-feature',
                    placeholder='Shared Key',
                    style={
                        'marginLeft': '2.5%',
                        'textAlign': 'center',
                        'width': '95%'
                    }
                ),

            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'marginLeft': '3%',
                    'width': '32%'
            }),
            html.Div([
                html.P('Merge-key: annot. table',
                    style={
                        'textAlign':'center',
                        'fontSize': 14
                    }),
                dcc.Dropdown(id='merge-key-annot',
                    placeholder='Shared Key',
                    style={
                        'marginLeft': '2.5%',
                        'textAlign': 'center',
                        'width': '95%'
                    }
                ),

            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'width': '32%'
            }),
            html.Div([
                html.P('Select external annots.',
                    style={
                        'textAlign':'center',
                        'fontSize': 14
                    }),
                dcc.Dropdown(id='external_annot',
                    placeholder='Annotations',
                    style={
                        'marginLeft': '2.5%',
                        'textAlign': 'center',
                        'width': '95%'
                    }
                ),


            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'marginTop': '0.5%',
                    'width': '32%'
            }),
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '340px',
                'width': '59%'
            }
        ),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.P('Generate and plot UMAP',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom':'0%'}),
        
        html.Div([
            html.P('UMAP Options',
                style={'textAlign': 'center',
                    'fontSize': 17,
                    'marginTop':'1.5%'}
            ),
            html.Hr(style={'marginTop':'0%', 'marginBottom':'2%'}),
            html.P('Feature scaling',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginTop':'0%',
                    'marginBottom': '0%'}
            ),
            dcc.Dropdown(id='feature-scaling',
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
                    'marginLeft': '5%',
                    'textAlign': 'center',
                    'width': '90%'
                }
            ),

            html.Hr(style={'marginTop':'4%', 'marginBottom':'4%'}),
            html.P('UMAP n_neighbors',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginBottom': '0%'}
            ),

            html.Div([
                dcc.Slider(
                    id='n_neighbors',
                    min=0,
                    max=100,
                    step=1,
                    value=10,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        25: '25',
                        50: '50',
                        100: '100',
                    },
                )],
                style={'width':'90%', 'marginLeft':'5%', 'marginTop':'3%'}
            ),
            html.Hr(style={'marginTop':'3%', 'marginBottom':'3%'}),

            html.P('UMAP min_dist',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginBottom': '0%'}
            ),

            html.Div([
                dcc.Slider(
                    id='min_dist',
                    min=0,
                    max=1,
                    step=0.05,
                    value=0.1,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        0.25: '0.25',
                        0.5: '0.5',
                        1: '1',
                    },
                )],
                style={'width':'90%', 'marginLeft':'5%', 'marginTop':'3%'}
            ),
            html.Hr(style={'marginTop':'3%', 'marginBottom':'3%'}),

            html.P('UMAP metric',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginBottom': '0%'}
            ),

            dcc.Input(
                    id = 'umap_metric',
                    type = 'text',
                    placeholder='euclidean',
                    style = {'width': '80%', 'marginTop': '1.5%',
                        'marginLeft': '10%'}
                ),

            html.Hr(style={'marginTop':'3%', 'marginBottom':'3%'}),

            html.P('K-Means Clustering',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginBottom': '0%'}
            ),
            html.P('Clustering is mutually exclusive with annotation color mapping.',
                style={'textAlign': 'center',
                    'fontSize': 12,
                    'marginTop':'0%',
                    'marginBottom': '0%',
                    'marginLeft':'10%',
                    'marginRight': '10%'}
            ),
            dcc.RadioItems(id='clustering',
                options=[
                    {'label': 'Cluster observations', 'value':'cluster'},
                    {'label': 'Do not cluster', 'value':'no_cluster'},

                ],
                value='no_cluster',
                style={
                    'marginLeft': '20%',
                    'textAlign': 'left',
                    'width': '90%'
                }
            ),

            html.Div([
                dcc.Slider(
                    id='n_cluster',
                    min=0,
                    max=30,
                    step=1,
                    value=10,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        10: '10',
                        20: '20',
                        30: '30',
                    },
                )],
                style={'width':'90%', 'marginLeft':'5%', 'marginTop':'3%'}
            ),
            html.P('# clusters',
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginTop': '1%'}
            ),
            html.Hr(style={'marginTop':'5%', 'marginBottom':'5%'}),
            html.Div([
                html.Button('Generate UMAP!', id = 'generate_umap',
                    style={
                    'width':'95%',
                    'marginLeft': '2.5%'
                    }
                )
            ],
                style={
                    'marginTop': '2%',
                    'marginLeft': '10%', 
                    'width': '80%',
                    'verticalAlign': 'top'
                }),
            html.Hr(style={'marginTop':'5%', 'marginBottom':'5%'}),

            dcc.Link("UMAP parameters docs",
                href="https://umap-learn.readthedocs.io/en/latest/parameters.html",
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginLeft':'10%'}
                ),
            html.Br(),
            dcc.Link("Feature scaling docs",
                href="https://scikit-learn.org/stable/auto_examples/preprocessing/plot_all_scaling.html",
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginLeft':'10%'}
                ),


        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '700px',
                'width': '20%',
                'borderRight': '1px #A6ACAF solid',
            }
        ),
        html.Hr()




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