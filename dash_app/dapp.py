import dash
from dash import html
import flask
from flask_caching import Cache
import pandas as pd
import requests


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


def query_panther(target_names, all_target_names, pval_thresh=0.1, enrichment_thresh=2, biological=True):
    url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep'
    # dataset ids from http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets


    panther_datasets = {
        "molecular_function": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
        "biological_process": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
        "cellular_component": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC"
    }
    if biological == 'bp':
        panth = panther_datasets['biological_process']
    elif biological == 'cc':

        panth = panther_datasets['cellular_component']
    else:
        panth = panther_datasets['molecular_function']

    panther_human_organism_id = 9606
    params = {
        'geneInputList': (','.join(target_names)),
        'refInputList': (','.join(all_target_names)),
        'organism': panther_human_organism_id,
        'refOrganism': panther_human_organism_id,
        'annotDataSet': panth,
        'enrichmentTestType': 'FISHER',
        'correction': 'FDR'
    }
    result = requests.post(url, params)
    d = result.json()
    # these are the hits (already sorted by p-value)
    annot = pd.DataFrame(data=d['results']['result'])
    annot['go_term_id'] = annot.term.apply(lambda s: s.get('id'))
    annot['go_term_label'] = annot.term.apply(lambda s: s.get('label'))
    annot.drop('term', axis=1, inplace=True)

    annot = annot[annot['pValue'] < pval_thresh]
    annot = annot[annot['fold_enrichment'] >= enrichment_thresh]
    annot.reset_index(drop=True, inplace=True)
    annot = annot[['number_in_list', 'expected', 'fold_enrichment', 'fdr',
        'pValue', 'go_term_id', 'go_term_label']].copy()

    return annot
