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
import plotly.graph_objects as go

fig = go.Figure()
fig.update_layout(template='simple_white')
fig.update_xaxes(showticklabels=False, ticks="")
fig.update_yaxes(showticklabels=False, ticks="")


def plotting_layout():
    return html.Div([
        html.Div([

            html.P('Use calculated enrichment table or upload a table',
                style={'textAlign': 'center',
                    'fontSize': 20,
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
                        'background-color': '#E5E7E9'}),

                html.P('Upload an enrichment or hits table', style={'textAlign': 'center',
                    'fontSize': 16, 'lineHeight': '15px', 'marginTop': '1%'}),
                dcc.Upload(
                    id='prep_table_upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select an enrichment or hits table')
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
                        'background-color': '#E5E7E9'},
                    # Only single file
                    multiple=False),


            ],
                style={
                    'display': 'inline-block',
                    'borderRight': '1px #A6ACAF solid',
                    'width': '49%',
                    'textAlign': 'center',
                    'vertical-align': 'top'}
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
                        ] ,

                        style={
                            'textAlign': 'center',
                            'width': '90%',
                            'marginLeft': '5%',
                            'marginTop': '3%'},

                        value='calculated'
                    ),
                ], style={
                    'display': 'inline-block',
                    'width': '65%',
                    'verticalAlign': 'top'}
                ),

                html.Div([
                    html.Button(
                        'Load data!',
                        id='load_button', style={'marginTop': '7%'})
                ], style={
                    'marginTop': '2%',
                    'display': 'inline-block',
                    'background-color': 'white',
                    'width': '35%'}),

                html.Hr(style={'marginTop': '2%', 'marginBottom': '3%'}),
                html.Div([
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
                            'background-color': '#E5E7E9'}),
                ], style={'display': 'inline-block',
                    'width': '50%', 'vertical-align': 'top'}),
                html.Div([
                    html.P('Select a marker label', style={'textAlign': 'center',
                        'fontSize': 16, 'lineHeight': '15px'}),

                    dcc.Dropdown(id='vol_marker_label',
                        placeholder='Select a label',
                        style={
                            'textAlign': 'center',
                            'marginTop': '2%',
                            'width': '90%',
                            'marginLeft': '5%'}
                    ),
                ], style={'display': 'inline-block',
                    'width': '50%', 'vertical-align': 'top'}),






            ],
                style={
                    'display': 'inline-block',
                    'width': '50%',
                    'vertical-align': 'top'}
            ),
        ]),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),

        html.Button('Call significant hits ▼',
            id='call_section', style={
                'border': '0px',
                'width': '80%',
                'marginLeft': '10%',
                'background-color': '#e8e8e8',
                'fontSize': 18}),


        html.Div([
            html.Div([
                html.P('Thresholding options', style={'textAlign': 'center',
                    'fontSize': 18, 'lineHeight': '15px', 'marginTop': '6%'}),

                dcc.RadioItems(id='thresh_option',
                    options=[
                        {'label': 'Hawaiian FDR (faster)',
                        'value': 'hawaii'},
                        {'label': 'Individual FDR',
                        'value': 'indiv'},
                        {'label': 'Manual FDR',
                        'value': 'manual'}
                    ],

                    style={
                        'textAlign': 'center',
                        'width': '90%',
                        'marginLeft': '5%',
                        'marginTop': '3%'
                    },
                    value='hawaii'
                ),
                html.P('Set FDR (%) threshold', style={'textAlign': 'center',
                    'fontSize': 18, 'lineHeight': '15px', 'marginTop': '4%'}),
                html.Div([
                    dcc.Slider(
                        id='fdr',
                        min=0,
                        max=40,
                        step=0.5,
                        value=5,
                        tooltip={"placement": "bottom", "always_visible": True},
                        marks={
                            0: '0',
                            10: '10',
                            20: '20',
                            30: '30',
                            40: '40'
                        },
                    )
                ],
                    style={'width': '90%', 'marginLeft': '5%', 'marginTop': '3%'}
                ),


            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'height': '200px',
                    'borderRight': '1px #A6ACAF solid',
                    'width': '24%'}
            ),

            html.Div([
                html.P('Set seeding threshold offset',
                    style={'textAlign': 'center',
                        'fontSize': 16,
                        'marginTop': '2%',
                        'marginBottom': '2%'}),

                html.Div([
                    dcc.Slider(
                        id='offset',
                        min=0,
                        max=15,
                        step=0.5,
                        value=2,
                        tooltip={"placement": "bottom", "always_visible": True},
                        marks={
                            0: '0',
                            5: '5',
                            10: '10',
                            15: '15'
                        },
                    )
                ],
                    style={'width': '95%', 'marginLeft': '2%', 'marginTop': '3%'}
                ),

                html.P('Set seeding threshold curvature',
                    style={'textAlign': 'center',
                        'fontSize': 16,
                        'marginTop': '2%',
                        'marginBottom': '2%'}),

                html.Div([
                    dcc.Slider(
                        id='curvature',
                        min=0,
                        max=10,
                        step=0.2,
                        value=3,
                        tooltip={"placement": "bottom", "always_visible": True},
                        marks={
                            0: '0',
                            5: '5',
                            10: '10',
                        },
                    )
                ],
                    style={'width': '95%', 'marginLeft': '3%', 'marginTop': '3%'}
                ),

            ],

                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'height': '200px',
                    'width': '27%'}
            ),

            html.Div([
                html.Button('Preview seeding FDR', id='preview',
                    style={
                        'fontSize': '13',
                        'marginLeft': '1%',
                        'width': '95%',
                        'marginTop': '8%',
                        'white-space': 'normal',
                        'background-color': 'white'}),

                html.Hr(style={'marginTop': '8%', 'marginBottom': '2%'}),
                html.P("estimated FDR:",
                    style={'textAlign': 'center',
                        'fontSize': 16,
                        'marginTop': '4%'}),
                html.P("  %", id='estimated_FDR',
                    style={'textAlign': 'center',
                        'fontSize': 16,
                        'marginTop': '2%'}),
            ], style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '200px',
                'width': '23%',
                'borderRight': '1px #A6ACAF solid'}
            ),




            html.Div([
                html.Button('Call significant hits!', id='hits_button',
                    style={
                        'marginLeft': '7.5%',
                        'width': '85%',
                        'marginTop': '15%',
                        'white-space': 'normal',
                        'background-color': 'white'}
                ),
                html.Button('Download hits table!', id='download_hits_button',
                    style={
                        'marginLeft': '7.5%',
                        'width': '85%',
                        'marginTop': '4%',
                        'white-space': 'normal'}
                ),


                dcc.Download(id="download_hits_table"),
            ],
                style={
                    'vertical-align': 'top',
                    'display': 'inline-block',
                    'height': '200px',
                    'width': '25%'}
            ),

        ], id='call_div', style={'display': 'none', 'marginTop': '1%'}),

        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),


        html.Div([

            html.P('Volcano plot',
                style={'textAlign': 'center',
                    'fontSize': 28}),
            html.Div([
                dcc.Dropdown(id='volcano_dropdown_1',
                    placeholder='Select a sample',
                    style={
                        'textAlign': 'center'
                    }
                )
            ], style={'display': 'inline-block', 'vertical-align': 'top',
                'width': '40%', 'marginLeft': '20%', 'marginBottom': '1%'}),
            html.Div([
                html.Button('Plot!', id='volcano_button_1',
                    style={
                        'marginLeft': '5%',
                        'vertical-align': 'top',
                        'width': '70%',
                        'white-space': 'normal'
                    }
                )
            ], style={'display': 'inline-block',
                'vertical-align': 'top',
                'width': '30%',
                'marginBottom': '1.5%'}),

            dcc.Graph(id='matrix_fig_1', figure=fig, style={'height': '90%'}),



        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '830px',
                'width': '53%'}),

        html.Div([

            dcc.Input(
                id='vol_search_plot',
                type='text',
                placeholder='Search label',
                style={'width': '100%', 'marginTop': '20%'}
            ),
            html.Button('Search!', id='vol_search_button',
                style={
                    'marginTop': '4%',
                    'width': '100%',
                    'white-space': 'normal'}),

            html.Hr(style={'marginTop': '5%', 'marginBottom': '5%',
                'lineWidth': 3}),

            dash_table.DataTable(
                id='vol_selected_data',
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
                    'height': '480px',
                    'width': '100%',
                    'border': '0.5px #BDC3C7 solid',
                },
                fixed_rows={'headers': True, 'data': 0}),
            html.P('0 points selected',
                id='vol_selection_count',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginTop': '15%'}),

            html.Button('Download', id='vol_download_subspace_button',
                style={
                    'marginTop': '4%',
                    'width': '100%',
                    'white-space': 'normal'}),



        ], style={
            'vertical-align': 'top',
            'display': 'inline-block',
            'height': '830px',
            'width': '12%'}),

        html.Div([
            html.P('Set manual FDR on a sample',
                style={'textAlign': 'center', 'fontSize': 20, 'marginTop': '0%'}),
            html.P('Offset',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginTop': '1%',
                    'marginBottom': '1%'}),

            html.Div([
                dcc.Slider(
                    id='man_offset',
                    min=0,
                    max=15,
                    step=0.5,
                    value=2,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        5: '5',
                        10: '10',
                        15: '15'
                    },
                )
            ],
                style={'width': '90%', 'marginLeft': '5%', 'marginTop': '1%'}
            ),

            html.P('Curvature',
                style={'textAlign': 'center',
                    'fontSize': 16,
                    'marginTop': '1%',
                    'marginBottom': '1%'}),

            html.Div([
                dcc.Slider(
                    id='man_curvature',
                    min=0,
                    max=10,
                    step=0.2,
                    value=3,
                    tooltip={"placement": "bottom", "always_visible": True},
                    marks={
                        0: '0',
                        5: '5',
                        10: '10',
                    },
                )
            ],
                style={'width': '90%', 'marginLeft': '5%', 'marginTop': '1%'}
            ),

            html.Button("Set FDR & plot", id='set_fdr_button',
                style={
                    'width': '60%',
                    'marginLeft': '20%',
                    'marginTop': '2%',
                }),

            html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),

            html.P('Plot options',
                style={'textAlign': 'center', 'fontSize': 20, 'marginTop': '1%'}),

            html.Div([
                dcc.Checklist(
                    id='plot_options',
                    options=[
                        {'label': 'Display FDR curve and interactors', 'value': 'fdr'},
                        {'label': 'Display labels for interactors', 'value': 'label'},
                        {'label': 'Highlight external annotations', 'value': 'ext'},
                    ],
                    value=[]),
            ], style={
                'marginLeft': '10%',
                'fontSize': 15,
                'marginTop': '3%'}),

            html.P('Select annotations',
                style={'textAlign': 'left',
                    'marginLeft': '13%',
                    'marginTop': '2%',
                    'fontSize': 18}),

            html.Div([

                dcc.Dropdown(id='vol_annot_select',
                    placeholder='Select annotations')

            ],
                style={
                    'textAlign': 'center',
                    'width': '80%',
                    'marginLeft': '10%',
                    'marginTop': '0.5%'}),

            html.Hr(style={'marginTop': '3%', 'marginBottom': '3%'}),

            html.P('External annotations', style={
                'marginTop': '1%',
                'fontSize': 20,
                'textAlign': 'center'}),

            dcc.Upload(
                id='vol_annot_table_upload',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select an annot. table')
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
                    'textAlign': 'center'}),

            html.Div([
                html.P('Merge-key: feature table',
                    style={
                        'textAlign': 'center',
                        'fontSize': 14}),

                dcc.Dropdown(id='vol_merge_key_feature',
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
                    'width': '47%'}),


            html.Div([
                html.P('Merge-key: annot. table',
                    style={
                        'textAlign': 'center',
                        'fontSize': 14}),

                dcc.Dropdown(id='vol_merge_key_annot',
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
                    'width': '48%'}),

            html.Div([
                html.P('Select external annots.',
                    style={
                        'textAlign': 'center',
                        'fontSize': 14}),

                dcc.Dropdown(id='vol_external_annot',
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
                    'marginTop': '2%',
                    'marginLeft': '25%',
                    'width': '50%'}),

            dcc.Input(
                id='vol_annot_label',
                placeholder='Annotation name',
                type='text',
                style={'width': '40%', 'marginTop': '2%', 'textAlign': 'center',
                    'display': 'inline-block', 'marginLeft': '10%'}
            ),

            html.Button("Link!", id='vol_merge_button',
                style={
                    'width': '35%',
                    'marginLeft': '5%',
                    'marginTop': '3%',
                    'display': 'inline-block'
                }),

        ],
            style={
                'vertical-align': 'top',
                'display': 'inline-block',
                'height': '830px',
                'marginLeft': '1%',
                'borderLeft': '1px #BDC3C7 solid',
                'width': '32%'}),

        html.Hr(style={'marginBottom': '2%', 'marginTop': '2%'}),

        html.Button('Annotation enrichment analysis ▼',
            id='vol_annot_section', style={
                'border': '0px',
                'width': '80%',
                'marginLeft': '10%',
                'background-color': '#e8e8e8',
                'fontSize': 18}),

        html.Div([
            html.Div([
                dcc.RadioItems(id='annot_analysis_opt',
                    options=[
                        {'label': 'Strip chart',
                        'value': 'strip'},
                        {'label': 'Box plot',
                        'value': 'box'}
                    ],
                    labelStyle={
                        'display': 'inline-block',
                        'marginRight': '3%'
                    },

                    style={
                        'textAlign': 'center',
                        'width': '45%',
                        'display': 'inline-block',
                        'marginLeft': '55%',
                        'marginTop': '1%'},

                    value='strip'
                )], style={'display': 'inline-block',
                        'vertical-align': 'top',
                        'width': '60%'}),
            html.Div([
                html.Button("Analyze!", id='strip_button',
                    style={
                        'width': '60%',
                        'marginLeft': '3%'})],
                style={
                'display': 'inline-block',
                'vertical-align': 'top',
                'width': '40%'}),


            html.Div([
                html.Div([
                    dcc.RangeSlider(
                        id='enrichment_range',
                        min=-20,
                        max=20,
                        step=0.5,
                        value=[-10, 10],
                        tooltip={"placement": "bottom", "always_visible": True},
                    )],
                    style={'width': '80%', 'marginLeft': '10%', 'marginTop': '3%'}
                ),

            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'width': '50%',
                'marginLeft': '10%',
                'marginTop': '2%'
            }
            ),
            html.Div([
                html.Button('Adjust range', id='range_button',
                    style={'textAlign': 'center',
                        'width': '60%',
                        'marginLeft': '3%'}),
            ], style={
                'display': 'inline-block',
                'verticalAlign': 'top',
                'marginTop': '2%',
                'width': '40%'
            }
            ),
            dcc.Graph(id='annot_fig', figure=fig),



        ], id='annotation_div', style={'display': 'none', 'marginTop': '2%'}),


        html.Hr(style={'marginBottom': '2%', 'marginTop': '2%'}),

        html.Button('Selection GO analysis ▼',
            id='vol_go_section', style={
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
                    dcc.Dropdown(id='vol_gene_selector',
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
                    dcc.Dropdown(id='vol_go_cat',
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
                    id='vol_pval_cutoff',
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
                    id='vol_enrich_cutoff',
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
                html.Button('GO analysis', id='vol_go_analysis',
                    style={
                        'width': '35%',
                        'marginLeft': '12.5%',
                        'marginRight': '5%',
                        'marginBottom': '2%'
                    }
                ),
                html.Button('Download GO table', id='vol_download_go_button',
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
                id='vol_go_top_table',
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

        ], id='vol_go_div', style={'display': 'none', 'marginTop': '2%'}),


        html.Hr(style={'marginTop': '2%', 'marginBottom': '2%'}),


        html.Div(
            id='hits_table', style={'display': 'none'}),

        dcc.Download(id="vol_download_subspace"),
        dcc.Download(id="vol_download_go"),



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
