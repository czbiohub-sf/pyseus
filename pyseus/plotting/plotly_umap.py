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
        matrix, node_name, cluster, opacity=0.8,
        width=800, height=600, highlight=None, unlabelled_color='#D0D3D4',
        unlabelled_opacity=0.1, hover_data=None):

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
            hover_data=hover_data,
            opacity=opacity,
            template='simple_white')
        fig.update_traces(marker=dict(size=5.5))

    else:
        labelling = matrix[cluster].isna()
        labelled = matrix[~labelling]
        unlabelled = matrix[labelling]
        unlabelled[cluster] = 'unlabelled'

        labelled.sort_values(by=cluster, inplace=True)

        fig1 = px.scatter(
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
            opacity=opacity,
            hover_data=hover_data,
            template='simple_white')
        fig1.update_traces(marker=dict(size=5.5))
        fig1.update(layout_coloraxis_showscale=False)


        fig2 = px.scatter(
            unlabelled,
            x='umap_1',
            y='umap_2',
            labels={
                'umap_1': 'UMAP 1',
                'umap_2': 'UMAP 2'
            },
            color=cluster,
            hover_name=node_name,
            opacity=unlabelled_opacity,
            hover_data=hover_data,
            color_discrete_sequence=[unlabelled_color],
            template='simple_white')

        fig = go.Figure(data=fig2.data + fig1.data)


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

    fig.update_xaxes(showticklabels=False, title_text='UMAP 1', ticks="")
    fig.update_yaxes(showticklabels=False, title_text='UMAP 2', ticks="")
    fig.update_layout(
        template='simple_white',
        legend=dict(
            font=dict(size=14)
        )
    )


    return fig
