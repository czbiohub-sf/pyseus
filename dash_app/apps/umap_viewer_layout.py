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
from dash.dash_table.Format import Format, Scheme, Trim


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
                    'marginLeft': '5%',
                    'width': '44%',
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
                    'width': '45%',
                    'vertical-align': 'top'}),

        ]),
        html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),


        html.Button('UMAP feature selection ▼',
            id='feature_section', style={
                'border': '0px',
                'width': '80%',
                'marginLeft': '10%',
                'background-color': '#e8e8e8',
                'fontSize': 18}),


        html.Div([
            html.Div([
                html.P('Loaded feature table dimensions',
                    style={'textAlign': 'center',
                        'fontSize': 20,
                        'marginBottom': '0%',
                        'marginTop': '13%'}),
                html.P('0 features X 0 observations, 0 annotations',
                    id='feature_dims',
                    style={'textAlign': 'center',
                        'fontSize': 18,
                        'marginBottom': '0%',
                        'marginTop': '1%'}),
                html.Div([
                    html.Button(
                        'Transpose matrix',
                        id='transpose_button',
                        style={'background-color': 'white'})
                ], style={
                    'marginTop': '2%',
                    'textAlign': 'center',
                    'width': '60%',
                    'marginLeft': '20%'}
                ),


            ], style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '290px',
                'marginLeft': '10%',
                'marginTop': '1%',
                'width': '35%',
                'borderRight': '1px #A6ACAF solid'}

            ),


            html.Div([
                html.P('Select features', style={'textAlign': 'center',
                    'fontSize': 20}),
                dcc.Checklist(
                    id='um_features_checklist',
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
                    'height': '290px',
                    'marginLeft': '3%',
                    'marginTop': '1%',
                    'width': '45%'}

            ),
        ], id='feature_div', style={'display': 'none', 'marginTop': '0.5%'}),


        html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),

        html.Button('UMAP Options ▼',
            id='umops_section', style={
                'border': '0px',
                'width': '80%',
                'marginLeft': '10%',
                'background-color': '#e8e8e8',
                'fontSize': 18}),

        html.Div([
            html.Div([
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
                    'marginLeft': '10%',
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
                    'width': '19%',
                    'marginLeft': '3%',
                    'marginRight': '3%',
                    'height': '100px'}

                ),
                html.Div([
                    html.P('UMAP metric',
                        style={'textAlign': 'center',
                            'fontSize': 15,
                            'marginBottom': '0%',
                            'marginTop': '5%'}),


                    dcc.Input(
                        id='umap_metric',
                        type='text',
                        value='euclidean',
                        style={'width': '80%', 'marginTop': '1%',
                            'marginLeft': '20%'}
                    ),

                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'borderLeft': '1px #A6ACAF solid',
                    'width': '18%',
                    'height': '100px'
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
                        style={'width': '80%', 'marginLeft': '10%', 'marginTop': '3%'}
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
                        style={'width': '80%', 'marginLeft': '10%', 'marginTop': '3%'}
                    ),

                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'width': '50%',
                }
                ),
            ], id='umops_div', style={'display': 'none', 'marginTop': '0.5%'}),


            html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

            html.Div([
                html.Button('Calculate UMAP!', id='generate_umap',
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
                    'marginLeft': '40%',
                    'marginBottom': '2%'}),

            dcc.Link("Feature scaling docs",
                href="https://scikit-learn.org/stable/auto_examples/preprocessing/plot_all_scaling.html",
                style={'textAlign': 'center',
                    'fontSize': 15,
                    'marginTop': '4%',
                    'marginLeft': '3%'}),

            html.Hr(style={'marginTop': '1%', 'marginBottom': '5%'}),

        ],
            style={
                'vertical-align': 'top',
                'width': '90%',
                'marginLeft': '5%'}

        ),
    ])


def plotting_layout():
    return html.Div([

        # upload component
        html.Div([
            html.P('UMAP table status',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom': '0%',
                    'marginTop': '1%'}),

            html.P('UMAP table not ready', id='umap_table_status',
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
            html.P('Upload a UMAP table', style={'textAlign': 'center',
                'fontSize': 20, 'lineHeight': '15px', 'marginTop': '1%'}),
            dcc.Upload(
                id='umap_table_upload',
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
                multiple=False),

            html.Button(
                'Upload table!',
                id='umap_load_button', style={'marginTop': '3%', 'width': '50%'})

        ], style={
            'display': 'inline-block',
            'height': '280px',
            'borderRight': '1px #A6ACAF solid',
            'marginLeft': '5%',
            'width': '45%',
            'textAlign': 'center',
            'vertical-align': 'top'}),

        html.Div([

            html.P('Link external annotations', style={
                'marginTop': '1%',
                'fontSize': 20,
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
            html.P('CSV/TSV files accepted. Table header MUST be in top-row.',
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

            ], style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'marginLeft': '3%',
                'width': '32%'}),


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
                        'width': '95%'}),

            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'width': '32%'}),

            html.Div([
                html.P('Select external annots.',
                    style={
                        'textAlign': 'center',
                        'fontSize': 14}),

                dcc.Dropdown(id='external_annot',
                    placeholder='Annotations',
                    style={
                        'marginLeft': '2.5%',
                        'textAlign': 'center',
                        'width': '95%'}),


            ], style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'marginTop': '0.5%',
                'width': '32%'}),


            dcc.Input(
                id='annot_label',
                placeholder='Annotation name',
                type='text',
                style={'width': '40%', 'marginTop': '2%', 'textAlign': 'center',
                    'display': 'inline-block', 'marginLeft': '10%'}
            ),


            html.Button("Link!", id='merge_button',
                style={
                    'width': '30%',
                    'marginLeft': '5%',
                    'marginTop': '2%',
                    'display': 'inline-block'
                }),
        ], style={
            'vertical-align': 'top',
            'display': 'inline-block',
            'height': '300px',
            'width': '47%'}),


        html.Button('Labels & Annotations ▼',
            id='label_section', style={
                'border': '0px',
                'width': '80%',
                'marginLeft': '10%',
                'background-color': '#e8e8e8',
                'fontSize': 18}),


        html.Div([
            html.Hr(style={'marginTop': '0%', 'marginBottom': '0%'}),


            html.Div([
                dcc.Tabs(
                    id='um_tabs',
                    value='data_select',
                    children=[
                        dcc.Tab(
                            label='Data Selection',
                            value='data_select',
                            children=html.Div([
                                html.P('x - axis',
                                    style={'textAlign': 'center',
                                        'fontSize': 16,
                                        'marginTop': '2%',
                                        'marginBottom': '1%'}),
                                html.Div([
                                    dcc.Dropdown(id='x_select',
                                        options=[
                                            {'label': 'UMAP 1', 'value': 'umap_1'},
                                            {'label': 'UMAP 2', 'value': 'umap_2'},
                                        ],
                                        value='umap_1')],
                                    style={
                                        'textAlign': 'center',
                                        'width': '80%',
                                        'marginLeft': '10%',
                                        'marginTop': '1%'}),
                                html.P('y - axis',
                                    style={'textAlign': 'center',
                                        'fontSize': 16,
                                        'marginTop': '2%',
                                        'marginBottom': '1%'}),
                                html.Div([
                                    dcc.Dropdown(id='y_select',
                                        options=[
                                            {'label': 'UMAP 1', 'value': 'umap_1'},
                                            {'label': 'UMAP 2', 'value': 'umap_2'},
                                        ],
                                        value='umap_2')],
                                    style={
                                        'textAlign': 'center',
                                        'width': '80%',
                                        'marginLeft': '10%',
                                        'marginTop': '1%'}),


                            ]),

                        ),
                        dcc.Tab(
                            label='Clustering',
                            value='clustering',
                            children=html.Div([

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
                                html.Button(
                                    'Cluster!',
                                    id='cluster_button', style={
                                        'marginTop': '3%',
                                        'width': '50%',
                                        'marginLeft': '25%'})
                            ]),
                        ),

                    ]
                ),

            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'height': '220px',
                    'marginLeft': '2%',
                    'width': '25%',
                    'borderRight': '1px #A6ACAF solid'}),

            html.Div([
                html.Div([
                    html.P('Select a marker label', style={'textAlign': 'center',
                            'fontSize': 20}),


                    html.Div([
                        dcc.Dropdown(id='um_label_select',
                            placeholder='Select a label')],
                        style={
                            'textAlign': 'center',
                            'width': '80%',
                            'marginLeft': '10%',
                            'marginTop': '1%'}),

                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'width': '50%'}),

                html.Div([
                    html.P('Select annotations',
                        style={'textAlign': 'center',
                            'fontSize': 20}),
                    html.Div([

                        dcc.Dropdown(id='annot_select',
                            placeholder='Select an annotation')],
                        style={
                            'textAlign': 'center',
                            'width': '80%',
                            'marginLeft': '10%',
                            'marginTop': '1%'}),
                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'width': '50%'}),

                html.Hr(style={'marginTop': '2.5%', 'marginBottom': '2.5%'}),

                html.P('Unlabelled marker style',
                    style={'textAlign': 'center',
                        'fontSize': 20,
                        'marginBottom': '0%'}),

                html.Div([
                    dcc.Link("Select color",
                        href="https://g.co/kgs/6JeS7i",
                        style={'textAlign': 'center',
                            'fontSize': 14,
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
                            'fontSize': 14,
                            'marginBottom': '0%',
                            'marginLeft': '7.5%',
                            'width': '90%'}),

                    dcc.Slider(
                        id='opacity',
                        min=0,
                        max=1,
                        step=0.05,
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
                ], style={
                    'display': 'inline-block',
                    'verticalAlign': 'top',
                    'width': '60%'
                }
                ),

            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'height': '220px',
                    'marginTop': '1%',
                    'width': '45%',
                    'borderRight': '1px #A6ACAF solid'}),

            html.Div([

                html.P('Download annotated table',
                    style={'textAlign': 'center',
                        'marginTop': '10%',
                        'fontSize': 20,
                        'marginBottom': '0%'}),
                dcc.Download(id='annot_table_dl'),
                html.Button(
                    'Download!',
                    id='annot_dl_button', style={'marginTop': '1%', 'width': '75%',
                        'marginLeft': '12.5%'}),

            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'height': '220px',
                    'marginTop': '1%',
                    'marginLeft': '1%',
                    'width': '20%'}),

        ], id='label_div', style={'display': 'none', 'marginTop': '2%'}),




        html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),


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
                'display': 'inline-block',
                'marginLeft': '2.5%',
                'height': '700px',
                'width': '80%'}),

        html.Div([
            html.Button(
                'Plot UMAP!',
                id='plot_button',
                style={'marginTop': '4%', 'width': '95%', 'marginLeft': '2.5%'}),

            html.Hr(style={'marginTop': '7%', 'marginBottom': '7%',
                'lineWidth': 3}),

            dcc.Input(
                id='search_plot',
                type='text',
                placeholder='Search label',
                style={'width': '95%', 'marginTop': '1.5%',
                    'marginLeft': '2.5%'}
            ),
            html.Button('Search!', id='search_button',
                style={
                    'marginTop': '4%',
                    'marginBottom': '4%',
                    'marginLeft': '2.5%',
                    'width': '95%',
                    'white-space': 'normal'}),

            html.Hr(style={'marginTop': '3%', 'marginBottom': '7%',
                'lineWidth': 3}),


            dash_table.DataTable(
                id='selected_data',
                columns=[{'name': 'Marker label', 'id': 'marker'}],
                page_size=10000,
                style_header={
                    'fontWeight': 'bold'
                },
                style_cell={
                    'textAlign': 'center'
                },
                style_table={
                    'overflowY': 'auto',
                    'overflowX': 'auto',
                    'height': '550px',
                    'width': '95%',
                    'border' : '0.5px #BDC3C7 solid'
                },
                fixed_rows={'headers': True, 'data': 0}),

            html.P('0 data points selected',
                id='selection_count',
                style={'textAlign': 'center',
                    'fontSize': 18,
                    'marginTop': '5%',
                    'marginBottom': '0.5%'}),
        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'marginLeft': '1%',
                'height': '700px',
                'width': '16%'}),



        html.Button('Download selection', id='download_subspace_button',
            style={
                'marginTop': '1%',
                'marginLeft': '30%',
                'width': '40%',
                'white-space': 'normal'}),
        html.Hr(style={'marginTop': '2%', 'marginBottom': '1%'}),

        html.Button('Selection GO analysis ▼',
            id='go_section', style={
                'border': '0px',
                'width': '80%',
                'marginLeft': '10%',
                'background-color': '#e8e8e8',
                'fontSize': 18}),

        html.Div([
            html.Div([
                html.P('Select a gene names column',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),
                html.Div([
                    dcc.Dropdown(id='gene_selector',
                        placeholder='Select gene names',
                        style={
                            'textAlign': 'center',
                            'width': '100%'}),
                ], style={
                    'marginLeft': '10%',
                    'width': '80%',
                    'marginTop': '1.5%'})

            ], style={
                'display': 'inline-block',
                'borderRight': '1px #A6ACAF solid',
                'marginLeft': '15%',
                'width': '20%',
                'textAlign': 'center',
                'vertical-align': 'top'}),

            html.Div([
                html.P('GO category',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),
                html.Div([
                    dcc.Dropdown(id='go_cat',
                        options=[
                            {'label': 'Cellular Component', 'value': 'cc'},
                            {'label': 'Biological Process', 'value': 'bp'},
                            {'label': 'Molecular Function', 'value': 'mf'}],

                        value='cc',
                        style={
                            'textAlign': 'center',
                            'width': '100%'}),
                ], style={
                    'marginLeft': '10%',
                    'width': '80%',
                    'marginTop': '1.5%'})

            ], style={
                'display': 'inline-block',
                'borderRight': '1px #A6ACAF solid',
                'width': '20%',
                'textAlign': 'center',
                'vertical-align': 'top'}),

            html.Div([
                html.P('GO p-val cutoff',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),

                dcc.Input(
                    id='pval_cutoff',
                    type='number',
                    value=0.1,
                    style={'width': '60%', 'marginTop': '1.5%', 'textAlign': 'center'}
                ),

            ], style={
                'display': 'inline-block',
                'borderRight': '1px #A6ACAF solid',
                'width': '15%',
                'textAlign': 'center',
                'vertical-align': 'top'}),

            html.Div([
                html.P('GO enrichment cutoff',
                    style={'textAlign': 'center',
                        'fontSize': 15,
                        'marginBottom': '0%'}),

                dcc.Input(
                    id='enrich_cutoff',
                    type='number',
                    value=2,
                    style={'width': '60%', 'marginTop': '1.5%', 'textAlign': 'center'}
                ),

            ], style={
                'display': 'inline-block',
                'width': '15%',
                'textAlign': 'center',
                'vertical-align': 'top'}),


            html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),


            html.Div([
                html.Button('GO analysis', id='go_analysis',
                    style={
                        'width': '35%',
                        'marginLeft': '12.5%',
                        'marginRight': '5%',
                        'marginBottom': '2%'
                    }
                ),
                html.Button('Download GO table', id='download_go_button',
                    style={
                        'width': '35%',
                        'white-space': 'normal'
                    }
                ),
            ],
                style={
                    'marginTop': '2%',
                    'marginLeft': '10%',
                    'width': '80%',
                    'verticalAlign': 'top'}),
            html.Hr(style={'marginTop': '1%', 'marginBottom': '1%'}),

            html.P('Top 10 enriched GO terms',
                style={'textAlign': 'center',
                    'fontSize': 20,
                    'marginBottom': '0%'}),

            dash_table.DataTable(
                id='go_top_table',
                columns=[
                    {'name': 'GO term', 'id': 'go_term'},
                    {'name': '# in selection', 'id': 'num'},
                    {'name': '# expected', 'id': 'expected'},
                    {'name': 'fold enrichment', 'id': 'fold'},
                    {'name': 'p-Val', 'id': 'pval', 'type': 'numeric',
                        'format': Format(precision=3, scheme=Scheme.exponent)},
                ],
                style_cell_conditional=[
                    {'if': {'column_id': 'num'},
                    'width': '15%'},
                    {'if': {'column_id': 'expected'},
                    'width': '15%'},
                    {'if': {'column_id': 'fold'},
                    'width': '15%'},
                    {'if': {'column_id': 'pval'},
                    'width': '15%'},
                    {'if': {'column_id': 'go_term'},
                    'width': '40%'},
                ],
                data=[{
                    'min': [0],
                    'max': [0],
                    'avg': [0],
                    'stdev': [0]
                }],
                style_header={
                    'fontWeight': 'bold'
                },
                style_cell={
                    'textAlign': 'center'
                },
                style_table={
                    'marginLeft': '7.5%',
                    'width': '85%',
                    'marginTop': '1%'
                },
                fixed_rows={'headers': True, 'data': 0}),

        ], id='go_div', style={'display': 'none', 'marginTop': '2%'}),


        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),




        dcc.Download(id="download_subspace"),
        dcc.Download(id="download_go"),







        # Hiddn divs inside the app for computations and storing intermediate values
        dcc.Store(id='um_cache_p_table'),


        html.Div(
            id='um_features', style={'display': 'none'}),
        html.Div(
            id='um_annots', style={'display': 'none'}),
        html.Div(
            id='final_features', style={'display': 'none'}),
        html.Div(
            id='table_dims', style={'display': 'none'}),
        html.Div(
            id='transposed_table', style={'display': 'none'}),
        html.Div(
            id='external_annot_series', style={'display': 'none'}),

    ])
