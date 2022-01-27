import dash
from dash import html
import flask 
from flask_caching import Cache


# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# server = flask.Flask(__name__)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
	suppress_callback_exceptions=True)

cache = Cache(app.server, config={
    'CACHE_TYPE': 'redis',
    # Note that filesystem cache doesn't work on systems with ephemeral
    # filesystems like Heroku.
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': 'cache-directory',

    # should be equal to maximum number of users on the app at a single time
    # higher numbers will store more data in the filesystem / redis cache
    'CACHE_THRESHOLD': 50
})




