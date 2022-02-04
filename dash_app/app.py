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
    'CACHE_THRESHOLD': 200
})


# def uploaded_processed_table(session_id, content=None):
#     """
#     Processed table cache for slots 1-3, and page-specific upload tables
#     """
#     @cache.memoize(args_to_ignore=['content'])
#     def upload_processed_table(session_id, content):
#         # parse file
#         content_type, content_string = content.split(',')
#         decoded = base64.b64decode(content_string)

#         raw_table = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
#             low_memory=False, header=[0,1], index_col=0)
        
#         return raw_table.to_json()
    
#     processed_table = upload_processed_table(session_id, content)
    
#     return pd.read_json(processed_table)
