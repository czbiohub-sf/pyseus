from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import uuid
import base64

from dapp import app
from dapp import saved_processed_table


import pandas as pd
import io

# import all pages in the app
from apps import preprocessing, matrix_viewer, umap_viewer, volcano_viewer, home
from apps import preprocessing_layout, matrix_viewer_layout, umap_viewer_layout,\
    volcano_calculation_layout, volcano_plotting_layout

session_id = str(uuid.uuid4())

server = app.server


# strings for initial population uploaded/saved table labels
no_table = 'No table uploaded'
no_p_table = 'No table processed'

first_link = {'fontSize': 20,
    'marginLeft': '14%',
    'text-decoration': 'none',
    'color': '#2F83C4',
    'fontWeight': 'bold'}
link_style = {'fontSize': 20,
    'marginLeft': '10%',
    'text-decoration': 'none',
    'color': '#2F83C4',
    'fontWeight': 'bold'}
section_style = {'width': '19%', 'textAlign': 'center', 'display': 'inline-block'}
first_style = {'marginLeft': '2.5%', 'width': '19%', 'textAlign': 'center',
    'display': 'inline-block'}

# Set up of navigation bar
navbar = html.Div([
    html.Div(style={'background-color': '#DCE7EC', 'height': '20px',
        'borderTop': '1px gray solid'}),

    html.Div([
        html.Div([
            dcc.Link("Home", href="/home", style=first_link),
        ], style=first_style),
        html.Div([
            dcc.Link("Pre-processing", href="/pre_processing", style=link_style),
        ], style=section_style),
        html.Div([
            dcc.Link("Clustergram", href="/clustergram", style=link_style),
        ], style=section_style),
        html.Div([
            dcc.Link("Volcano Plot", href="/volcano", style=link_style),
        ], style=section_style),
        html.Div([
            dcc.Link("UMAP", href="/umap", style=link_style),
        ], style=section_style),

    ], style={'background-color': '#DCE7EC',
        'height': '50px',
        'borderBottom': '1px gray solid'})
])

index_layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content'),

    # Slot labels for uploaded / saved tables that are cached server-side
    # Acts as 'global' labels as they are saved in the index layout
    html.Div(children=no_table, id='slot_label_1', style={'display': 'none'}),
    html.Div(children=no_table, id='slot_label_2', style={'display': 'none'}),
    html.Div(children=no_table, id='slot_label_3', style={'display': 'none'}),
    html.Div(children=no_p_table, id='slot_label_4', style={'display': 'none'}),
    html.Div(children=no_p_table, id='slot_label_5', style={'display': 'none'}),
    html.Div(children=no_p_table, id='slot_label_6', style={'display': 'none'}),

    dcc.Store(data=session_id, id='session_id')

])

# embedding the navigation bar
app.layout = index_layout

# all the layouts for validation layout
pp_upload_layout = preprocessing_layout.upload_layout()
pp_process_layout = preprocessing_layout.process_layout()
matrix_calc_layout = matrix_viewer_layout.calculation_layout()
matrix_plotting_layout = matrix_viewer_layout.plotting_layout()
umap_custom_layout = umap_viewer_layout.customize_layout()
umap_plotting_layout = umap_viewer_layout.plotting_layout()
volc_calc_layout = volcano_calculation_layout.calculation_layout()
volc_plot_layout = volcano_plotting_layout.plotting_layout()

# "complete" layout for validation
app.validation_layout = html.Div([
    index_layout,
    pp_upload_layout,
    pp_process_layout,
    matrix_calc_layout,
    matrix_plotting_layout,
    umap_custom_layout,
    umap_plotting_layout,
    volc_calc_layout,
    volc_plot_layout,
    home.layout
])

# open and cache preload ground truths
itzhak = pd.read_csv('preload_annots/dan_organelle_truths_20220707.csv')
manu = pd.read_csv('preload_annots/MLgroup_organelle_curation_2.4.csv')

saved_processed_table(session_id + 'itzhak', itzhak)
saved_processed_table(session_id + 'manu', manu)




@app.callback(Output('page-content', 'children'),
        Input('url', 'pathname'))
def display_page(pathname):
    if pathname == '/pre_processing':
        return preprocessing.layout
    elif pathname == '/clustergram':
        return matrix_viewer.layout
    elif pathname == '/volcano':
        return volcano_viewer.layout
    elif pathname == '/umap':
        return umap_viewer.layout
    else:
        return home.layout


if __name__ == "__main__":
    app.run_server(debug=True)
