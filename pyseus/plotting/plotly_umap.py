import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Normalizer


from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import plotly.offline
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import plotly.express as px
import umap

def scale_table(matrix, method):
    """
    takes a feature table and scale the data accordingly
    """


    if method == 'standard':
        scaler = StandardScaler()
    elif method == 'robust':
        scaler = RobustScaler()
    elif method == 'min_max':
        scaler = MinMaxScaler()
    elif method == 'l1_norm':
        scaler = Normalizer(norm='l1')
    elif method == 'l2_norm':
        scaler = Normalizer(norm='l2')
    else:
        # if there is no specified scaler, return values with nulls dropped
        return matrix

    scaled = scaler.fit_transform(matrix)

    return scaled


def interaction_umap(
        matrix, node_name, cluster, opacity=0.5,
        width=800, height=600, highlight=None):

    matrix = matrix.copy()
    matrix.reset_index(inplace=True, drop=False)
    
    if node_name == 'None':
        node_name = None

    if cluster == 'None':
        fig = px.scatter(
            matrix,
            x='umap_1',
            y='umap_2',
            labels={
            'umap_1': 'UMAP 1',
            'umap_2': 'UMAP 2'
            },
            hover_name=node_name,
            opacity=opacity)
        fig.update_traces(marker=dict(size=5.5))       

    else: 
        labelling = matrix[cluster].isna()
        labelled = matrix[~labelling]
        unlabelled = matrix[labelling]

        labelled.sort_values(by=cluster, inplace=True)

        fig = px.scatter(
            labelled,
            x='umap_1',
            y='umap_2',
            labels={
            'umap_1': 'UMAP 1',
            'umap_2': 'UMAP 2'
            },
            hover_name=node_name,
            color=cluster,
            color_continuous_scale=px.colors.cyclical.mygbm[: -1],
            opacity=opacity)
        fig.update_traces(marker=dict(size=5.5))
        fig.update(layout_coloraxis_showscale=False)
        fig.add_scatter(
            x=unlabelled['umap_1'],
            y=unlabelled['umap_2'],
            mode='markers',
            showlegend=False,
            hoverinfo='skip',
            opacity=0.2,
            marker=dict(color='grey'))

        # if highlight:
        #     labelled = matrix[matrix[node_name].isin(highlight)]
        #     fig.add_scatter(
        #         x=labelled['umap_1'],
        #         y=labelled['umap_2'],
        #         mode='markers',
        #         showlegend=False,
        #         hoverinfo='text',
        #         opacity=1,
        #         text=labelled[node_name],
        #         marker=dict(color='#fc4f30', size=14))

    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(
        legend=dict(
            font=dict(size=14)
        )
        )
    
    return fig


