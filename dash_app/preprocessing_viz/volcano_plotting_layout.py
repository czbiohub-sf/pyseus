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