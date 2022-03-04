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


def customize_layout():
    return html.Div([


        html.Div([

            html.P('Upload a processed feature table or use a pre-loaded one',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom': '0%'}),
            html.P("""The UMAP generating script requires a standard format
                of feature table. Please use the Process Table page to create one.""",
                style={
                    'textAlign': 'center',
                    'fontSize': 14,
                    'marginTop': '0%'}),

            # upload component
            html.Div([
                dcc.Upload(
                    id='um_raw_table_upload',
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
                # html.P('Uploaded feature table', id='um_raw_table_filename',
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
                    html.Button('Read table!', id='um_read_table_button')],
                    style={
                        'marginTop': '2%',
                        'marginLeft': '30%',
                        'width': '40%',
                        'verticalAlign': 'top'}
                ),
            ],
                style={
                    'display': 'inline-block',
                    'borderRight': '1px #A6ACAF solid',
                    'width': '49%',
                    'textAlign': 'center',
                    'vertical-align': 'top'}
            ),
            html.Div([

                dcc.Dropdown(id='um_preloaded_dropdown',
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
                        id='um_preload_button',)
                ], style={
                    'marginTop': '2%',
                    'textAlign': 'center',
                    'display': 'inline-block',
                    'width': '40%',
                    'marginLeft': '30%'

                }),

            ],
                style={
                    'display': 'inline-block',
                    'width': '50%',
                    'vertical-align': 'top'}),

        ]),
        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
        html.P('Loaded feature table dimensions',
            style={'textAlign': 'center',
                'fontSize': 20,
                'marginBottom': '0%',
                'marginTop': '1%'}),
        html.P('0 features X 0 observations, 0 annotations',
            id='feature_dims',
            style={'textAlign': 'center',
                'fontSize': 16,
                'marginBottom': '0%',
                'marginTop': '0%'}),
        html.Div([
            html.Button(
                'Transpose matrix',
                id='transpose_button',
                style={'background-color': 'white'})
        ], style={
            'marginTop': '1%',
            'textAlign': 'center',
            'width': '40%',
            'marginLeft': '30%'}
        ),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),


        html.Div([
            html.Div([
                html.P('Select a marker label', style={'textAlign': 'center',
                    'fontSize': 18}),
            ],
                style={
                    'vertical-align': 'top',
                    'marginLeft': '10%',
                    'width': '80%'}

            ),
            dcc.Dropdown(id='um_label_select',
                placeholder='Select a label',
                style={
                    'textAlign': 'left',
                    'width': '90%',
                    'marginLeft': '5%'
                }
            ),
            html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),

            html.P('Select features for UMAP generation', style={'textAlign': 'center',
                'fontSize': 18}),
            dcc.Checklist(
                id='um_features_checklist',
                style={
                    'overflowY': 'auto',
                    'overflowX': 'auto',
                    'height': '210px',
                    'marginLeft': '10%',
                    'width': '80%',
                    'border': '0.5px #BDC3C7 solid',
                }),
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '360px',
                'marginLeft': '3%',
                'marginTop': '1%',
                'width': '24%',
                'borderRight': '1px #A6ACAF solid'}

        ),
        html.Div([
            html.Div([
                html.P('Annotation options for UMAP color-mapping',
                    style={'textAlign': 'center',
                        'fontSize': 18}),
            ],
                style={
                    'vertical-align': 'top',
                    'marginLeft': '10%',
                    'width': '80%'}

            ),

            dcc.RadioItems(id='annot_options',
                options=[
                    {'label': 'No annotations', 'value': 'no_annot'},
                    {'label': 'Internal annotations', 'value': 'internal'},
                    {'label': 'External annotations', 'value': 'external'},
                    {'label': 'Clustering (K-Means)', 'value': 'cluster'},

                ],
                value='no_annot',
                style={
                    'marginLeft': '5%',
                    'textAlign': 'left',
                    'width': '90%'
                }
            ),

            html.Hr(style={'marginTop': '6%', 'marginBottom': '6%'}),

            html.P('K-Means clustering',
                style={'textAlign': 'center',
                    'fontSize': 18,
                    'marginBottom': '0%'}),


            html.P('Designate # of clusters',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginTop': '2%',
                    'marginBottom': '2%'}),

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
                style={'width': '90%', 'marginLeft': '5%', 'marginTop': '3%'}
            ),



        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '360px',
                'marginTop': '1%',
                'width': '25%',
                'borderRight': '1px #A6ACAF solid'}),


        html.Div([

            html.P('Internal annotations',
                style={'textAlign': 'center',
                    'fontSize': 18,
                    'marginTop': '2%'}),
            html.Div([
                dcc.Dropdown(id='annot_select',
                    placeholder='Select a label',
                    style={
                        'textAlign': 'center'
                    }
                ),
            ], style={
                'vertical-align': 'top',
                'marginLeft': '20%',
                'width': '60%'}),

            html.Hr(style={'marginTop': '3%', 'marginBottom': '0%'}),
            html.P('External annotations', style={
                'marginTop': '1%',
                'fontSize': 18,
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
                    'background-color': '#E5E7E9'},

                # Only single file
                multiple=False
            ),
            html.P('csv or tsv files accepted. Table header MUST be in top-row.',
                style={
                    'fontSize': 13,
                    'textAlign': 'center'
                }),
            html.Div([
                html.P('Merge-key: feature table',
                    style={
                        'textAlign': 'center',
                        'fontSize': 14
                    }),
                dcc.Dropdown(id='merge_key_feature',
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
                        'textAlign': 'center',
                        'fontSize': 14
                    }),
                dcc.Dropdown(id='merge_key_annot',
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
                        'textAlign': 'center',
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
            html.Button("Link external annotations", id='merge_button',
                style={
                    'width': '60%',
                    'marginLeft': '20%',
                    'marginTop': '3%'
                }),
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '360px',
                'width': '47%'}

        ),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),
    ])


def plotting_layout():
    return html.Div([

        html.Div([
            html.P('UMAP Options',
                style={'textAlign': 'center',
                    'fontSize': 24,
                    'marginTop': '1.5%'}),


            html.Hr(style={'marginTop': '0%', 'marginBottom': '1%'}),
            html.Div([
                html.P('Feature scaling',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginTop': '5%',
                        'marginBottom': '0%'}),

                dcc.Dropdown(id='feature_scaling',
                    options=[
                        {'label': 'None', 'value': 'None'},
                        {'label': 'Standardization', 'value': 'standard'},
                        {'label': 'Robust Scaler', 'value': 'robust'},
                        {'label': 'Min-max scaling', 'value': 'min_max'},
                        {'label': 'l1 normalization', 'value': 'l1_norm'},
                        {'label': 'l2 normalization', 'value': 'l2_norm'},
                    ],
                    value='standard',
                    style={
                        'marginLeft': '5%',
                        'textAlign': 'center',
                        'width': '90%'
                    }
                ),
            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'width': '24%',
                'borderRight': '1px #A6ACAF solid',
                'height': '100px'}

            ),
            html.Div([
                html.P('Random State',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),

                html.P('Enter an integer',
                    style={'textAlign': 'center',
                        'fontSize': 12,
                        'marginBottom': '0%',
                        'marginLeft': '7.5%',
                        'width': '85%'}),


                dcc.Input(
                    id='random_state',
                    type='text',
                    value='None',
                    style={'width': '80%', 'marginTop': '1.5%',
                        'marginLeft': '10%'}
                ),
            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'width': '18%',
                'borderRight': '1px #A6ACAF solid',
                'height': '100px'}

            ),
            html.Div([
                html.P('UMAP metric',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%',
                        'marginTop': '3%'}),


                dcc.Input(
                    id='umap_metric',
                    type='text',
                    value='euclidean',
                    style={'width': '80%', 'marginTop': '1%',
                        'marginLeft': '10%'}
                ),

            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'width': '18%',
                'borderRight': '1px #A6ACAF solid',
                'height': '100px'
            }
            ),
            html.Div([
                html.P('Unlabelled markers (for annots.)',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),

                html.Div([
                    dcc.Link("Select color",
                        href="https://g.co/kgs/6JeS7i",
                        style={'textAlign': 'center',
                            'fontSize': 13,
                            'marginBottom': '0%',
                            'marginLeft': '40%',
                            'width': '30'}),


                    dcc.Input(
                            id='marker_color',
                            type='text',
                            value='#D0D3D4',
                            style={'width': '80%', 'marginTop': '1.5%',
                                'marginLeft': '10%'}
                    ),
                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'width': '40%'
                }
                ),
                html.Div([
                    html.P('Opacity',
                        style={'textAlign': 'center',
                            'fontSize': 13,
                            'marginBottom': '0%',
                            'marginLeft': '7.5%',
                            'width': '90%'}),

                    dcc.Slider(
                        id='opacity',
                        min=0,
                        max=1,
                        step=0.05,
                        value=0.1,
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
                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'width': '60%'
                }
                ),


            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'height': '100px',
                'width': '38%'
            }
            ),

            html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),
            html.Div([
                html.P('UMAP n_neighbors',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),


                html.Div([
                    dcc.Slider(
                        id='n_neighbors',
                        min=0,
                        max=20,
                        step=1,
                        value=5,
                        tooltip={"placement": "bottom", "always_visible": True},
                        marks={
                            0: '0',
                            5: '5',
                            10: '10',
                            15: '15',
                            20: '20'
                        },
                    )],
                    style={'width': '90%', 'marginLeft': '5%', 'marginTop': '3%'}
                ),

            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'width': '49%',
                'borderRight': '1px #A6ACAF solid',
            }
            ),

            html.Div([

                html.P('UMAP min_dist',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),


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
                    style={'width': '90%', 'marginLeft': '5%', 'marginTop': '3%'}
                ),

            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'width': '50%',
            }
            ),


            html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),

            html.Div([
                html.Button('Generate UMAP!', id='generate_umap',
                    style={
                        'width': '35%',
                        'marginLeft': '10%',
                        'marginRight': '10%',
                        'marginBottom': '2%'
                    }
                ),
                html.Button('Download UMAP', id='download_umap_button',
                    style={
                        'width': '35%',
                        'white-space': 'normal'
                    }
                ),
                dcc.Download(id="download_umap"),
            ],
                style={
                    'marginTop': '2%',
                    'marginLeft': '10%',
                    'width': '80%',
                    'verticalAlign': 'top'}),

            dcc.Link("UMAP parameters docs",
                href="https://umap-learn.readthedocs.io/en/latest/parameters.html",
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginTop': '4%',
                    'marginLeft': '25%',
                    'marginBottom': '2%'}),

            dcc.Link("Feature scaling docs",
                href="https://scikit-learn.org/stable/auto_examples/preprocessing/plot_all_scaling.html",
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginTop': '4%',
                    'marginLeft': '3%'}),

            html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),

        ],
            style={
                'vertical-align': 'top',
                'width': '90%',
                'marginLeft': '5%'}

        ),

        html.P('UMAP Figure',
            style={'textAlign': 'center',
                'fontSize': 24,
                'marginTop': '1.5%'}),


        html.Div([
            dcc.Graph(id='umap_fig',
                style={'height': '100%'})
        ],
            style={
                'vertical-align': 'top',
                'marginLeft': '2.5%',
                'height': '700px',
                'width': '95%'}),

        html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),

        html.P('Select UMAP subspace for download',
            style={'textAlign': 'center',
                'fontSize': 24,
                'marginTop': '1.5%',
                'marginBottom': '1.5%'}),

        html.Div([
            html.P('x-range',
                style={'textAlign': 'center',
                    'fontSize': 18,
                    'marginBottom': '0%'}),


            html.Div([
                dcc.RangeSlider(
                    id='umap_x',
                    min=0,
                    max=20,
                    step=1,
                    value=[5, 15],
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        5: '5',
                        10: '10',
                        15: '15',
                        20: '20'
                    },
                )],
                style={'width': '80%', 'marginLeft': '10%', 'marginTop': '3%'}
            ),

        ], style={
            'display': 'inline-block',
            'verticalAlign': 'top',
            'width': '49%',
            'borderRight': '1px #A6ACAF solid'}),


        html.Div([
            html.P('y-range',
                style={'textAlign': 'center',
                    'fontSize': 18,
                    'marginBottom': '0%'}),


            html.Div([
                dcc.RangeSlider(
                    id='umap_y',
                    min=0,
                    max=20,
                    step=1,
                    value=[5, 15],
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        5: '5',
                        10: '10',
                        15: '15',
                        20: '20'
                    },
                )],
                style={'width': '80%', 'marginLeft': '10%', 'marginTop': '3%'}
            ),

        ], style={
            'display': 'inline-block',
            'verticalAlign': 'top',
            'width': '50%'}),

        html.Button('Download Subspace', id='download_subspace_button',
            style={
                'marginTop': '3%',
                'marginLeft': '30%',
                'width': '40%',
                'white-space': 'normal'
            }
        ),
        dcc.Download(id="download_subspace"),

        html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),





        # Hiddn divs inside the app for computations and storing intermediate values
        dcc.Store(id='um_cache_p_table'),

        html.Div(
            id='um_table', style={'display': 'none'}),
        html.Div(
            id='um_processed_table', style={'display': 'none'}),
        html.Div(
            id='um_features', style={'display': 'none'}),
        html.Div(
            id='um_annots', style={'display': 'none'}),
        html.Div(
            id='table_dims', style={'display': 'none'}),
        html.Div(
            id='transposed_table', style={'display': 'none'}),
        html.Div(
            id='external_annot_series', style={'display': 'none'}),

    ])
