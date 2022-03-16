import dash
from dash import html
import flask
from flask_caching import Cache
import pandas as pd


# initiate app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

server = flask.Flask(__name__)

app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True)
app.enable_dev_tools(debug=True, dev_tools_hot_reload=False)

# set up server-side cache.

cache = Cache(app.server, config={
    "CACHE_DEFAULT_TIMEOUT": 43200,
    # Note that filesystem cache doesn't work on systems with ephemeral
    # filesystems like Heroku.
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': 'cache-directory',

    # should be equal to maximum number of users on the app at a single time
    # higher numbers will store more data in the filesystem / redis cache
    'CACHE_THRESHOLD': 200
})


def saved_processed_table(session_id, processed_table=None, overwrite=False):
    """
    save data tables to server-side cache using unique IDs found in different app pages

    """


    @cache.memoize(args_to_ignore=['processed_table'])
    def save_processed_table(session_id, processed_table):

        return processed_table.to_json()

    if overwrite:
        # delete cache to save a new one if overwrite option is True
        try:
            cache.delete_memoized(save_processed_table, session_id)

        except AttributeError:
            pass

    processed_table = save_processed_table(session_id, processed_table)

    return pd.read_json(processed_table)


def cycle_style_colors(style, color_1='#DCE7EC', color_2='#dcdfec'):
    """
    cycle between two colors of button
    """
    if style is None:
        style = {}

    if 'background-color' not in style.keys():
        style['background-color'] = color_1

    elif style['background-color'] == color_1:
        style['background-color'] = color_2
    else:
        style['background-color'] = color_1

    return style
