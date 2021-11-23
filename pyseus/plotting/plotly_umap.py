import pandas as pd
import requests
import numpy as np
import initial_processing as ip

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import plotly.offline
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import plotly.express as px
import umap

def org_balance_scale_split_predict(table, ref_col, quant_cols, max_sample=150):
    table = table.copy()
    all_cols = quant_cols + [ref_col]
    table = table[all_cols]

    # find the lowest value count of a label
    table.dropna(inplace=True)
    num_labels = table[ref_col].nunique()
    refs = table[ref_col].unique()
    refs.sort()

    samples = []
    for col in refs:
        sample = table[table[ref_col]==col]
        if sample.shape[0] <= max_sample:
            samples.append(sample)
        else:
            random_sample = sample.sample(max_sample)
            samples.append(random_sample)
    
    balanced = pd.concat(samples).reset_index(drop=True)

    labels = pd.factorize(balanced[ref_col])
    definitions = labels[1]

    balanced['label'] = labels[0]

    X = balanced[quant_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    y = balanced['label'].values

    # split and standard scale
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2)


    # Random Forest Classifier
    classifier = RandomForestClassifier()
    classifier.fit(X_train, y_train)
    
    y_pred = classifier.predict(X_test)
     
    # Reverse factorizers
    reversefactor = dict(zip(range(num_labels), definitions))
    
    y_tested = np.vectorize(reversefactor.get)(y_test)
    y_predicted = np.vectorize(reversefactor.get)(y_pred)

    return y_tested, y_predicted



def one_balance_scale_split_predict(table, ref_col, quant_cols, ref):
    table = table.copy()
    all_cols = quant_cols + [ref_col]
    table = table[all_cols]

    # find the lowest value count of a label
    table.dropna(inplace=True)
    refs = table[ref_col].unique()
    refs.sort()


    sample = table[table[ref_col]==ref]
    sample_size = sample.shape[0]
    others = table[table[ref_col]!=ref]
    # balancing
    others = others.sample(sample_size)
    others[ref_col] = 'Others'

    
    balanced = pd.concat([sample, others]).reset_index(drop=True)

    labels = pd.factorize(balanced[ref_col])
    definitions = labels[1]

    balanced['label'] = labels[0]

    X = balanced[quant_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    y = balanced['label'].values

    # split and standard scale
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2)


    # Random Forest Classifier
    classifier = RandomForestClassifier()
    classifier.fit(X_train, y_train)
    
    y_pred = classifier.predict(X_test)
     
    # Reverse factorizers
    reversefactor = dict(zip(range(2), definitions))
    
    y_tested = np.vectorize(reversefactor.get)(y_test)
    y_predicted = np.vectorize(reversefactor.get)(y_pred)

    return y_tested, y_predicted, classifier


def confusion_precision_recall(tests, predictions, exp=''):
    cross_recall = pd.crosstab(
        tests, predictions, rownames=['Actual Compartment'],
            colnames=['Predicted Compartment']).apply(
                lambda r: np.round(r/r.sum(),3), axis=1)
    
    cross_precision = pd.crosstab(
        tests, predictions, rownames=['Actual Compartment'],
            colnames=['Predicted Compartment']).apply(
                lambda r: np.round(r/r.sum(),3), axis=0)
    
    orgs = list(cross_recall)
    recall_table = {}
    precision_table = {}
    for org in orgs:
        recall_table[org] = cross_recall[org][org]
        precision_table[org] = cross_precision[org][org]
    
    recall_table = pd.DataFrame(recall_table, index=[exp])
    precision_table = pd.DataFrame(precision_table, index=[exp])

    return cross_recall, recall_table, precision_table



def interaction_umap(
        matrix, quant_cols, node_name, cluster, n_neighbors=3, min_dist=0.5, metric='euclidean', opacity=0.5,
        width=800, height=600, highlight=None):

    matrix = matrix.copy()
    matrix.dropna(subset=quant_cols,inplace=True)
    X = matrix[quant_cols].values

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    fit = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric
    )

    u = fit.fit_transform(X_scaled)
    matrix['umap_1'] = u[: , 0]
    matrix['umap_2'] = u[: , 1]
    matrix.reset_index(inplace=True, drop=False)

    labelling = matrix[cluster].isna()
    labelled = matrix[~labelling]
    unlabelled = matrix[labelling]

    fig = px.scatter(
        labelled,
        x='umap_1',
        y='umap_2',
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
    if highlight:
        labelled = matrix[matrix[node_name].isin(highlight)]
        fig.add_scatter(
            x=labelled['umap_1'],
            y=labelled['umap_2'],
            mode='markers',
            showlegend=False,
            hoverinfo='text',
            opacity=1,
            text=labelled[node_name],
            marker=dict(color='#fc4f30', size=14))

    fig.update_layout(
        width=width,
        height=height,
        margin={'l': 30, 'r': 30, 'b': 20})
    fig.show()


