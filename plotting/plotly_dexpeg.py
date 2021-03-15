import matplotlib.pyplot as plt
import matplotlib
from numbers import Number
import numpy as np
import pandas as pd
import plotly.offline
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import plotly.express as px
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler



def spatial_kmeans_cluster(spatial_df, n_clusters=8):
    spatial_df = spatial_df.copy()

    fractions = ['01K', '03K', '05K', '12K', '24K', '80K', 'Cyt',
        'Ratio_005', 'Ratio_Nuclease', 'Ratio_No_Nuclease',
        'Upper_norm', 'Lower_norm']

    sdata = spatial_df[fractions].values
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)

    cluster_labels = KMeans(n_clusters=n_clusters).fit_predict(scaled)
    col_label = 'cluster_' + str(n_clusters)

    spatial_df[col_label] = cluster_labels

    return spatial_df


def dexpeg_umap(spatial_df, label='', opacity=0.5, width=800, height=600,
        show_ticks=False, misc_hover=False):
    """
    Ways of showing interactive HEK spatial maps
    """
    spatial_df = spatial_df.copy()
    # fractions = ['1K', 'Cytosol', 'Input', 'Upper_norm', 'Lower_norm']
    # for fraction in fractions:
    #     spatial_df[fraction] = spatial_df[fraction].apply(np.round, args=[3])

    if 'search' in label:
        search = label.split('_')[1].lower()
        spatial_df['Gene names'] = spatial_df['Gene names'].astype(str)
        labelling = spatial_df['Gene names'].map(lambda x: search in x.lower())
        labelled = spatial_df[labelling]
        unlabelled = spatial_df[~labelling]
        fig = px.scatter(
            labelled,
            x='umap_1',
            y='umap_2',
            hover_name='Gene names',
            # hover_data=fractions,
            opacity=opacity)
        fig.update_traces(marker=dict(size=5.5, color='#008fd5'))
        if misc_hover:
            hover = 'text'
        else:
            hover = 'skip'

        fig.add_scatter(
            x=unlabelled['umap_1'],
            y=unlabelled['umap_2'],
            mode='markers',
            showlegend=False,
            hoverinfo=hover,
            opacity=0.2,
            text=unlabelled['Gene names'],
            marker=dict(color='grey'))
    elif label == 'manu_label':
        # mode for showing reference genes
        # divide the into labeled and unlabelled
        labelling = spatial_df[label].apply(lambda x: True if
            type(x) == str else False)
        labelled = spatial_df[labelling]

        unlabelled = spatial_df[~labelling]
        labelled.sort_values(by=label, inplace=True)
        unlabelled[label] = unlabelled[label].apply(
            lambda x: 'Unlabelled')
        fig = px.scatter(
            labelled,
            x='umap_1',
            y='umap_2',
            hover_name='Gene names',
            # hover_data=['Nuc', '03K', '05K', '12K', '24K', '80K', 'Cyt'],
            color=label,
            color_discrete_sequence=px.colors.qualitative.Light24,
            # color_discrete_sequence=px.colors.cyclical.mrybm,
            opacity=opacity)
        fig.update_traces(marker=dict(size=6))
        if misc_hover:
            hover = 'text'
        else:
            hover = 'skip'

        fig.add_scatter(
            x=unlabelled['umap_1'],
            y=unlabelled['umap_2'],
            mode='markers',
            showlegend=False,
            hoverinfo='skip',
            opacity=0.2,
            text=unlabelled['Gene names'],
            marker=dict(color='grey'))
    elif 'cluster' in label:
        spatial_df[label] = spatial_df[label].apply(str)
        spatial_df.sort_values(by=label, inplace=True)
        fig = px.scatter(
            spatial_df,
            x='umap_1',
            y='umap_2',
            color=label,
            color_discrete_sequence=px.colors.qualitative.T10,
            hover_name='Gene names',
            # hover_data=['Nuc', '03K', '05K', '12K', '24K', '80K', 'Cyt'],
            opacity=opacity)

    elif 'feature' in label:
        label = label.split('_', 1)[1]
        fig = px.scatter(
            spatial_df,
            x='umap_1',
            y='umap_2',
            color=label,
            color_continuous_scale='Geyser',
            hover_name='Gene names',
            # hover_data=['Nuc', '03K', '05K', '12K', '24K', '80K', 'Cyt'],
            opacity=opacity)

    elif label == 'reference':
        # mode for showing reference genes
        # divide the into labeled and unlabelled
        labelling = spatial_df[label].apply(lambda x: True if
            type(x) == str else False)
        labelled = spatial_df[labelling]

        unlabelled = spatial_df[~labelling]
        labelled.sort_values(by=label, inplace=True)
        unlabelled[label] = unlabelled[label].apply(
            lambda x: 'Unlabelled')
        fig = px.scatter(
            labelled,
            x='umap_1',
            y='umap_2',
            hover_name='Gene names',
            # hover_data=['Nuc', '03K', '05K', '12K', '24K', '80K', 'Cyt'],
            color=label,
            color_discrete_sequence=px.colors.qualitative.Light24,
            # color_discrete_sequence=px.colors.cyclical.mrybm,
            opacity=opacity)
        fig.update_traces(marker=dict(size=6))
        if misc_hover:
            hover = 'text'
        else:
            hover = 'skip'

        fig.add_scatter(
            x=unlabelled['umap_1'],
            y=unlabelled['umap_2'],
            mode='markers',
            showlegend=False,
            hoverinfo='skip',
            opacity=0.2,
            text=unlabelled['Gene names'],
            marker=dict(color='grey'))

    else:
        fig = px.scatter(
            spatial_df,
            x='umap_1',
            y='umap_2',
            hover_name='Gene names',
            color='NOC_prediction',
            # hover_data=fractions,
            opacity=opacity)
        if misc_hover:
            hover = 'text'
        else:
            hover = 'skip'

    fig.update_xaxes(showticklabels=show_ticks, title_text='')
    fig.update_yaxes(showticklabels=show_ticks, title_text='')
    fig.update_layout(
        width=width,
        height=height,
        title={'text': "Integrated Spatial Proteomics",
        'x': 0.5},
        margin={'l': 30, 'r': 30, 'b': 20})
    fig.show()
    # fig.add_scatter(
    #     x=unlabelled['umap_1'],
    #     y=unlabelled['umap_2'],
    #     mode='markers',
    #     showlegend=False,
    #     hoverinfo='skip',
    #     opacity=0.2,
    #     text=unlabelled['Gene names'],
    #     marker=dict(color='grey'))


def get_top_corrs(df, gene_name, top_n=20, protein_id='', metric='corr', ascending=False):
    df = df.copy()
    
    gene_name = gene_name.upper()
    # check if gene name is in the list
    if gene_name not in list(df):
        print(gene_name + " not in the map!")
        return

    # check if multiple gene names in the selection
    selection = df[df['Gene names'] == gene_name]
    if selection.shape[0] > 1:
        if len(protein_id) == 0:
            print('Multiple protein groups with identical Gene names')
            print('Please specify a protein ID')
            print('Available Protein IDs:')
            print(selection['Protein IDs'].to_list())
        else:
            selection = selection[selection['Protein IDs'] == protein_id]
            if selection.shape[0] == 0:
                print("Invalid Protein ID")
                return
    if 'Protein IDs' in list(selection):
        selection = selection.drop(columns=['Protein IDs', 'Gene names']).T
    else:
        selection = selection.drop(columns=['Gene names']).T
    selection.columns = [metric]
    selection = selection.sort_values(by=metric, ascending=ascending).iloc[:top_n]
    return selection
