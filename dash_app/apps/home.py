import base64
import datetime
import markdown
import io
import json
import simplejson
import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dcc
from dash import html
from dash import dash_table


layout = html.Div([
    html.P('Welcome to the Pyseus Explorer',
        style={'textAlign': 'center', 'fontSize': 28, 'marginTop':'2%',
            'marginBottom': '1%'}),
    html.Hr(),
    html.P('This web-app package includes data exploratory tools specific for mass-spec outputs.\
        Please use the drop-down menu on the top-right to navigate to different tools.\
        For documentation on how to use the tools, please visit here.',
        style={'textAlign': 'center', 'fontSize': 18, 'marginTop':'2%',
            'marginBottom': '1%', 'marginLeft': '15%', 'width':'70%'}),
    html.P('Developed by Kibeom Kim, for internal use at CZ Biohub only.',
        style={'textAlign': 'center', 'fontSize': 16, 'marginTop':'2%',
            'marginBottom': '1%', 'marginLeft': '15%', 'width':'70%'}),

])