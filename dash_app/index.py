from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from app import app
# import all pages in the app
from apps import preprocessing, matrix_viewer, umap_viewer, volcano_viewer, home




navbar = html.Div([
html.Div(style={'background-color': '#DCE7EC', 'height':'20px', 'borderTop': '1px gray solid'}),
html.Div([
    html.Div([
        dcc.Link("Home", href="/home", style={'fontSize': 20, 'marginLeft':'10%'}),
    ], style={'marginLeft': '10%', 'width':'16%', 'textAlign':'top', 'display': 'inline-block'}),

    html.Div([
        dcc.Link("Pre-processing", href="/pre_processing", style={'fontSize': 20, 'marginLeft':'10%'}),
    ], style={'width':'16%', 'textAlign':'top', 'display': 'inline-block'}),


   html.Div([
        dcc.Link("Clustergram", href="/clustergram", style={'fontSize': 20, 'marginLeft':'10%'}),
    ], style={'width':'16%', 'textAlign':'top', 'display': 'inline-block'}),


    html.Div([
        dcc.Link("Volcano Plot", href="/volcano", style={'fontSize': 20, 'marginLeft':'10%'}),

    ], style={'width':'16%', 'textAlign':'top', 'display': 'inline-block'}),
    html.Div([
        dcc.Link("UMAP", href="/umap", style={'fontSize': 20, 'marginLeft':'10%'}),
    ], style={'width':'16%', 'textAlign':'top', 'display': 'inline-block'})

],
style={'background-color': '#DCE7EC', 'height':'50px', 'borderBottom': '1px gray solid'})
])
# embedding the navigation bar
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content')
])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
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