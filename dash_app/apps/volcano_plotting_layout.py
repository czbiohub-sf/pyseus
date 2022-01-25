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

def plotting_layout():
    return html.Div([
        html.Div([
            
            html.P('Use calculated enrichment table or upload a table',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom':'0%',
                    'marginTop': '3%',
                    'marginBottom': '1%'}),

            # upload component
            html.Div([
                html.P('Enrichment table not calculated', id='enrichment_table_status', 
                    style={
                        'marginBottom': '3%',
                        'marginLeft': '10%',
                        'width': '80%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'solid',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'background-color': '#E5E7E9'
                    }
                ),
            html.P('Upload an enrichment or hits table', style={'textAlign': 'center',
                'fontSize': 16, 'lineHeight':'15px', 'marginTop':'1%'}),
                dcc.Upload(
                    id='prep_table_upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select a table')
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

                
                # html.Div([
                #     html.Button('Read enrichment table!', id='read_enrichment_button'),

                #     html.Button('Read hits table!', id='read_hits_button',
                #         style={
                #             'marginLeft':'4%'
                #         })],
                #         style={
                #             'marginTop': '2%',
                #             'marginLeft': '5%', 
                #             'width': '90%',
                #             'verticalAlign': 'top',
                #             'white-space': 'normal'
                #         }),
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
                html.Div([
                    dcc.RadioItems(id='load_options',
                        options=[
                            {'label': 'Use calculated enrichment table',
                            'value': 'calculated'},
                            {'label': 'Use uploaded enrichment table',
                            'value': 'pre_enrichment'},
                            {'label': 'Use uploaded hits table',
                                'value': 'pre_hits'}
                            ],
                
                        style={
                            'textAlign': 'center',
                            'width': '90%',
                            'marginLeft': '5%',
                            'marginTop':'3%'
                        },
                        value='calculated'
                    ),
                ], style={
                    'display': 'inline-block',
                    'width': '65%',
                    'verticalAlign':'top'
                }
                ),

                html.Div([
                    html.Button(
                        'Load data!',
                        id =  'load_button', style={'marginTop': '7%'})
                    ], style={
                        'marginTop': '2%',
                        'display': 'inline-block',

                        'width': '35%'

                    }),

                html.Hr(style={'marginTop': '2%', 'marginBottom': '3%'}),
                html.P('Select a marker label', style={'textAlign': 'center',
                    'fontSize': 16, 'lineHeight':'15px'}),
                dcc.Dropdown(id='vol_marker_label',
                    placeholder='Select a label',
                    style={
                        'textAlign': 'center',
                        'marginTop': '2%',
                        'width': '70%',
                        'marginLeft': '15%',
                    }
                ),



            ],
            style = {
                'display': 'inline-block',
                'width': '50%',
                'vertical-align': 'top',
                }
            ),
        ]),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.P('Call significant hits',
                        style={'textAlign': 'center',
                            'fontSize': 20,
                            'marginTop':'1%'}
                ),
        html.Div([
           
            html.P('Thresholding options', style={'textAlign': 'center',
                'fontSize': 18, 'lineHeight':'15px', 'marginTop':'6%'}),

            dcc.RadioItems(id='thresh_option',
                options=[
                    {'label': 'Hawaiian FDR (faster)',
                    'value': 'hawaii'},
                    {'label': 'Individual FDR',
                    'value': 'indiv'}
                    ],
        
                style={
                    'textAlign': 'center',
                    'width': '90%',
                    'marginLeft': '5%',
                    'marginTop': '3%'
                },
                value='indiv'
            ),
            html.P('Set FDR (%) threshold', style={'textAlign': 'center',
                'fontSize': 18, 'lineHeight':'15px', 'marginTop':'8%'}),
            html.Div([
                dcc.Slider(
                    id='fdr',
                    min=0,
                    max=10,
                    step=0.5,
                    value=5,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        5: '5',
                        10: '10',
                    },
                )],
                style={'width':'90%', 'marginLeft':'5%', 'marginTop':'3%'}
            ),

            
            ],

            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '200px',
                'borderRight': '1px #A6ACAF solid',
                'width': '24%',
            }
        ),

        html.Div([
            html.P('Set seeding threshold offset',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginTop': '2%',
                    'marginBottom': '2%'}
            ),
            html.Div([
                dcc.Slider(
                    id='offset',
                    min=0,
                    max=10,
                    step=0.5,
                    value=2,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        5: '5',
                        10: '10',
                    },
                )],
                style={'width':'90%', 'marginLeft':'5%', 'marginTop':'3%'}
            ),

            html.P('Set seeding threshold curvature',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginTop': '2%',
                    'marginBottom': '2%'}
            ),
            html.Div([
                dcc.Slider(
                    id='curvature',
                    min=0,
                    max=5,
                    step=0.2,
                    value=3,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        2: '2',
                        4: '4',
                    },
                )],
                style={'width':'90%', 'marginLeft':'5%', 'marginTop':'3%'}
                ),

            ],

            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '200px',
                'width': '40%',

                'borderRight': '1px #A6ACAF solid',
            }
        ),

        html.Div([
            html.Button('Call significant hits!', id='hits_button',
                style={
                    'marginLeft':'7.5%',
                    'width': '85%',
                    'marginTop': '4%',
                    'white-space': 'normal'
                }
            ),
            html.Button('Download hits table!', id='download_hits_button',
                style={
                    'marginLeft':'7.5%',
                    'width': '85%',
                    'marginTop': '4%',
                    'white-space': 'normal'
                }
            ),
            html.P('Hits table not ready', id='hits_table_status', 
                style={
                    'marginTop': '4%',
                    'marginLeft': '7.5%',
                    'width': '85%',
                    'height': '40px',
                    'lineHeight': '40px',
                    'borderWidth': '1px',
                    'borderStyle': 'solid',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'background-color': '#E5E7E9'
                }
            ),
            dcc.Download(id="download_hits_table"),
        ],
        style={
            'vertical-align': 'top',
            'display': 'inline-block',
            'height': '200px',
            'width': '35%',

        }
        ),
        
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),


        html.Div([
            html.P('Volcano plot (1)',
                style={'textAlign': 'center',
                    'fontSize': 18}
            ),
            dcc.Dropdown(id='volcano_dropdown_1',
                placeholder='Select a sample',
                style={
                    'textAlign': 'center',
                    'width': '90%',
                    'marginLeft': '5%',
                }
            ),
            html.Button('Plot volcano!', id='volcano_button_1',
                style={
                    'marginLeft':'20%',
                    'width': '60%',
                    'marginTop': '2%',
                    'white-space': 'normal'
                }
            ),
            dcc.Graph(id='matrix_fig_1', style=
                {'height': '100%'})
 
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '700px',
                'width': '49%',
                'borderRight': '1px #A6ACAF solid',
            }
        ),
        html.Div([
            html.P('Volcano plot (2)',
                style={'textAlign': 'center',
                    'fontSize': 18}
            ),
            dcc.Dropdown(id='volcano_dropdown_2',
                placeholder='Select a sample',
                style={
                    'textAlign': 'center',
                    'width': '90%',
                    'marginLeft': '5%',
                }
            ),
            html.Button('Plot volcano!', id='volcano_button_2',
                style={
                    'marginLeft':'20%',
                    'width': '60%',
                    'marginTop': '2%',
                    'white-space': 'normal'
                }
            ),
            dcc.Graph(id='matrix_fig_2', style=
            {'height': '100%'})
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '700px',
                'width': '50%',
            }),


        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.Div(
            id='hits_table', style={'display': 'none'})





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