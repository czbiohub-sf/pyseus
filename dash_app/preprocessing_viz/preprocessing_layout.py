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
            style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'2%',
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
                        html.A('Select a MS table')
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
                html.P('Uploaded MS table', id='raw_table_filename', 
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
                    dcc.RadioItems(
                        id='raw_table_sep',
                        options=[
                            {'label': 'comma-separated', 
                            'value': ','},
                            {'label': 'tab-separated',
                            'value': '\t'}
                        ],
                        value='\t',
                        ),
                    ],
                    style={
                        'display': 'inline-block',
                        'marginTop': '0.5%',
                        'width': '30%'
                    }),
                html.Div([
                    html.P(
                        '# top rows to skip', 
                        style={
                            'fontSize': 14,
                            'textAlign' : 'center',
                            'marginTop' : '0.5%',
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
                        'marginTop': '0.5%',
                        'width': '23%'
                    }),
                html.Div([
                    html.Button('Read table!', id = 'read_table_button')],
                    style={
                        'display': 'inline-block',
                        'marginTop': '2%',
                        'width': '30%',
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
                dcc.Upload(
                    id='config_upload',
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
                html.P('Uploaded configs file', id='config_filename', 
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
                    html.Button(
                        'Pre-load configs',
                        id =  'load-config-button',)
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


        # Column selection
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),    
        html.P('Select sample and metadata columns',
                style={'textAlign': 'center', 'fontSize': 22}),
        
        html.Div([
            html.Div([
                html.P('Sample columns', style={'textAlign': 'center',
                    'fontSize': 16}),
                dcc.Checklist(
                    id='sample_cols_checklist',
                    style={
                        'marginLeft': '10%',
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '200px',
                        'width': '80%',
                        'border': '0.5px #BDC3C7 solid',
                    }),
                html.P('Regular expression search', style = {'textAlign': 'center',
                    'marginTop': '1.5%'}),
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
            html.P(id='num_samples_selected',
                style={'textAlign': 'center',
                    'marginTop': '2%'})

            ],
            style = {
                'display': 'inline-block',
                'width': '49%',
                'borderRight': '1px #A6ACAF solid'
                }
            ),
            html.Div([
                html.P('Metadata columns', style={'textAlign': 'center',
                    'fontSize':16}),
                dcc.Checklist(
                    id='meta_cols_checklist',
                    style={
                        'marginLeft': '10%',
                        'overflowY': 'auto',
                        'overflowX': 'auto',
                        'height': '270px',
                        'width': '80%',
                        'border': '0.5px #BDC3C7 solid',
                    }),
                html.P(id='num_meta_selected',
                    style={'textAlign': 'center', 'marginTop': '3%'})],
                style = {
                    'display': 'inline-block',
                    'width': '50%',
                    'vertical-align': 'top'}
            ),    
        ]), 


        html.Div([
            html.Button(
                'Save sample / metadata columns',
                id=  'save-cols-button',
                style = {
                    'marginTop': '3%',
                    'width': '40%'
                }),
            ], style={
                'textAlign': 'center'
            }
        ),

        # Processing options
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.P('MS table processing options',
            style={'textAlign': 'center', 'fontSize': 22, 'marginBottom': '2%'}),     
        html.Div([
            dcc.Checklist(
                id='filter_rows_check',
                options=[{'label': 'Filter rows (MaxQuant-specific)',
                    'value': 'filter_maxquant'}],
                style={'marginTop': '5%'}
            ),
            html.P("""Filter potential contaminants, reverse-sequences,
                and protein groups only identified by site.""", style={
                    'fontSize': 13,
                    'width': '95%',
                    'marginBottom': '3%'
                }),
            dcc.Checklist(
                id='rename_samples_check',
                options=[{'label': 'Rename samples',
                    'value': 'rename_samples'}]),
            html.P("""IMPORTANT for standard sample naming format, which is
                used in downstream analyses and app functions. 
                Please use the renaming preview at the bottom of the page
                to specify regular expressions""", style={
                    'fontSize': 13,
                    'width': '95%'
                }),


        ], style={
            'display': 'inline-block',
            'marginLeft': '10%',
            'width': '39%',
            'borderRight': '1px #A6ACAF solid',
            'height': '220px',
            'vertical-align': 'top'

        }),   
       html.Div([
            html.Div([
                dcc.Checklist(
                    id='log_transform_check',
                    options=[{'label': 'Log transform data',
                        'value': 'log_transform'}]),


            ], style={'display': 'inline-block', 'width': '40%',
                    'vertical-align': 'top', 'marginTop': '3%'}),
            html.Div([
                dcc.RadioItems(
                    id='transform_option',
                    options=[
                        {'label': 'Log 2', 'value': 'log2'},
                        {'label': 'Log 10', 'value': 'log10'}
                    ],
                    value='log2',
                    labelStyle={
                        'display': 'inline-block',
                        'fontSize': 13,
                        'marginRight': '2%'
                    }
                ),
                
            ], style={'display': 'inline-block', 'width': '50%',
                'vertical-align': 'top', 'marginTop': '3%'}),
            html.P("Values of 0 will be substituted as nulls.", style={
                'fontSize': 13
            }),

            dcc.Checklist(
                id='remove_incomplete_rows',
                options=[{'label': 'Remove incomplete rows',
                    'value': 'remove_rows'}],
                style={'marginTop': '3%'}),
            html.P("""Remove rows that do not have at least one
                sample group that has real values in every technical replicate.""",
                style={
                    'fontSize': 13,
                    'width': '95%',
                    'marginBottom': '3%'
                }),
            html.Div([
                dcc.Checklist(
                    id='merge_reps',
                    options=[{'label': 'Merge technical replicates',
                        'value': 'merge_reps'}]),


            ], style={'display': 'inline-block', 'width': '70%',
                    'vertical-align': 'top'}),
            html.Div([
                dcc.RadioItems(
                    id='merge_option',
                    options=[
                        {'label': 'mean', 'value': 'mean'},
                        {'label': 'median', 'value': 'median'}
                    ],
                    value= 'median',
                    labelStyle={
                        'display': 'inline-block',
                        'fontSize': 14,
                        'marginRight': '2%'
                    }
                ),
                

            ], style={'display': 'inline-block', 'width': '30%',
                'vertical-align': 'top'}),
            html.P("Calculate average statistic among replicates. Null values are ignored.", style={
                'fontSize': 13,
                'marginBottom': '3%'
            }),

        ], style={
            'display': 'inline-block',
            'marginLeft': '4%',
            'width': '40%',
            'height': '220px',
            'vertical-align': 'top'
        }),   
        dcc.Checklist(
            id='replace_nulls',
            options=[{'label': 'Replace null values',
                'value': 'replace_nulls'}],
            style={
                'marginTop': '2%',
                'textAlign': 'center'
            }
           
            ),
            
        dcc.RadioItems(
            id='replace_option',
            options=[
                {'label': 'sub in zeroes', 'value': 'zero'},
                {'label': 'sample-wise imputation', 'value': 'sample_impute'},
                {'label': 'row-wise imputation', 'value': 'row_impute'}
            ],
            value='sample_impute',

            style={
                'textAlign': 'center',
                'marginTop': '0.7%'
            },

            labelStyle={
                'display': 'inline-block',
                'fontSize': 13,
                'marginRight': '1.5%'
            }
            ),
        html.P('Imputation distance',
            style={
                'textAlign': 'center',
                'fontSize': 14,
                'marginTop': '1%'
            }),
        html.Div([
            dcc.Slider(
                id='imputation_dist',
                min=0,
                max=3.5,
                step=0.1,
                value=1.8,
                tooltip={"placement": "bottom", "always_visible": True},
                marks={
                    0: '0',
                    0.5: '0.5',
                    1: '1',
                    1.5: '1.5',
                    2: '2',
                    2.5: '2.5',
                    3: '3',
                    3.5: '3.5'
                },
            )
        ], style={'width': '60%', 'marginLeft':'20%'}),

        html.P('Imputation width',
            style={
                'textAlign': 'center',
                'fontSize': 14,
                'marginTop': '1%'
            }),
        html.Div([
            dcc.Slider(
                id='imputation_width',
                min=0,
                max=1,
                step=0.1,
                value=0.3,
                tooltip={"placement": "bottom", "always_visible": True},
                marks={
                    0: '0',
                    0.2: '0.2',
                    0.4: '0.4',
                    0.6: '0.6',
                    0.8: '0.8',
                    1: '1'
                },
            )
        ], style={'width': '60%', 'marginLeft':'20%'}),
        html.P("""Imputation distance and width are in units of standard deviations
            of either the sample or row distributions.""",
            style={
                'textAlign': 'center',
                'fontSize': 13,
                'marginTop': '1%',
                'marginLeft': '20%',
                'width': '60%'
            }),

        # Process and Save
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.P('Run MS table processing and save files',
            style={'textAlign': 'center', 'fontSize': 22, 'marginBottom': '1.5%'}), 
        html.Button(
            'Process MS table!',
            id=  'process_table_button',
            style = {
                'marginLeft': '30%',
                'width': '40%',
                'height': '4%'
            }),
        html.P('Please make sure to double check all options.',
            style={'textAlign': 'center', 'fontSize': 13, 'marginTop':'0.5%'}),

        html.Div([
            dcc.Input(
                    id = 'ms_save_name',
                    type = 'text',
                    placeholder='ms table savefile name',
                    style = {'width': '80%', 'marginTop': '1.5%'}
                ),
            html.P('ex: 20211105_CZBSU01_preprocessed_table.csv',
                style={
                    'fontSize': 13,
                    'marginLeft': '2%'
                }),
            html.Button('Download processed MS table',
                id='download_table',
                style={
                    'textAlign': 'center',
                    'width': '70%',
                    'marginTop': '3%'
                })

        ],
        style={
            'display': 'inline-block',
            'textAlign': 'center',
            'marginLeft': '10%',
            'width': '39%',
            'borderRight': '1px #A6ACAF solid',
            'height': '10%',
            'vertical-align': 'top'

        }),

       html.Div([
            dcc.Input(
                    id = 'configs_save_name',
                    type = 'text',
                    placeholder='configs savefile name',
                    style = {'width': '80%', 'marginTop': '1.5%'}
                ),
            html.P('ex: 20211105_CZBSU_configs.csv',
                style={
                    'fontSize': 13,
                    'marginLeft': '2%'
                }),
            html.Button('Download configs file',
                id='download_configs',
                style={
                    'textAlign': 'center',
                    'width': '70%',
                    'marginTop': '3%'
                })
        ],
        style={
            'display': 'inline-block',
            'textAlign': 'center',
            'width': '40%',
            'height': '10%',
            'vertical-align': 'top'
        }),


        
        # Preview Renaming
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.P('Renaming samples preview',
            style={'textAlign': 'center', 'fontSize': 22}),
        html.P("""Sample names should follow the format of 'experiment_sample_rep#'""",
            style={'textAlign':'center', 'marginBottom':'0%'}),
        html.P("""ex: CZBSU01_LAMP1_1 or 20211031_infected-24hr_2""", 
            style={'textAlign':'center'}),
        html.Div([
            html.Div([
                dash_table.DataTable(
                    id='original_sample_cols',
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
                        'width': '80%',
                        'border' : '0.5px #BDC3C7 solid'
                    },
                    fixed_rows={'headers':True, 'data':0})],
                style = {
                    'display': 'inline-block',
                    'width': '50%',
                    'height': '250px'}
            ),
            html.Div([
                dash_table.DataTable(
                    id='renamed_sample_cols',
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
                        'width': '80%',
                        'border' : '0.5px #BDC3C7 solid'
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
                'width': '35%'}),
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
                'width': '55%'}),

        ]),
        html.Div([
            html.Div([
                html.P("Replacement regular expressions, split by semicolons (type NONE for deletion):")   
            ], style={
                'display': 'inline-block',
                'marginLeft': '5%',
                'width': '35%'}),
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
                'width': '55%'}),

        ]),
        html.Button('Test renaming samples', id = 'preview_rename',
            style={
                'width': '30%',
                'marginLeft': '35%',
                'marginTop': '2%'

            }),
        html.Div([
            html.Span(
                id='re_warning',
                style={'textAlign': 'center',
                'color':'red'})
            ]),
        html.Hr(),
            

        # Hiddn divs inside the app for computations and storing intermediate values
        html.Div(
            id='raw_table', style={'display': 'none'}),
        html.Div(
            id='processed_table', style={'display': 'none'}),
        html.Div(
            id='all_cols', style={'display': 'none'}),
        html.Div(
            id='sample_cols', style={'display': 'none'}),
        html.Div(
            id='meta_cols', style={'display': 'none'}),
        dcc.Download(id="download-ms-table-csv"),
        dcc.Download(id="download-configs-csv"),
        # html.Div(
        #     id='renamed_sample_cols', children=None, style={'display': 'none'}),
        # html.Div(
        #     id='renamed_table', children=None, style={'display': 'none'}),
        # html.Div(
        #     id='grouped_table', children=None, style={'display': 'none'}),
        
    ])
