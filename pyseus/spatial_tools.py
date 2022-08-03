import multiprocessing
import sys
import itertools
import scipy
import random
import re
import pandas as pd
import numpy as np
import anndata as ad


from pyseus import basic_processing as bp
from pyseus import primary_analysis as pa
from pyseus import validation_analysis as va
from external import clustering_workflows

from multiprocessing import Queue
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from itertools import repeat
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform
from scipy.stats import percentileofscore
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier



class SpatialTables():
    """
    SpatialTables class takes as input standard_table from
    AnalysisTables class and performs various processing tools
    that enhance the datasets
    """

    def __init__(self, hit_table, target_col='target', prey_col='prey',
            control_mat=None):
        """
        initiate class with a standard hit table
        """
        self.hit_table = hit_table
        self.target = target_col
        self.prey = prey_col
        self.control_mat = control_mat
        self.create_enrichment_table()

    def create_enrichment_table(self):
        """
        many of the spatial functions require an enrichment table,
        which is a wide-version of the hits table with just enrichments (excluding pvals)
        """

        # create an AnalysisTables class which handles the enrichment conversion
        analysis = pa.AnalysisTables()
        analysis.simple_pval_table = self.hit_table
        analysis.convert_to_enrichment_table(enrichment='enrichment', simple_analysis=True)

        self.enrichment_table = analysis.enrichment_table.copy()

    def enrichment_corr_control_mat(self, corr=0.4):
        """
        create a control matrix that has a correlation filter between samples,
        effectively removing samples that are closely related for pval-calculations
        """

        enrichments = self.enrichment_table.copy()
        enrichments = enrichments['sample'].copy()
        if self.control_mat is None:
            print("Please assign control_mat to the Class!")
            return
        else:
            control_mat = self.control_mat
        cols = list(enrichments)

        # get correlation tables
        sample_corrs = enrichments[cols].corr()

        # get the list of columns in mat table
        mat_samples = control_mat['Samples'].to_list()

        # create a new mat and filter contrasts by correlation
        new_mat = control_mat.copy()
        for col in cols:
            truths = []
            passing_corrs = sample_corrs[col][sample_corrs[col] <= corr]
            for sample in mat_samples:
                if sample in passing_corrs:
                    truths.append(True)
                else:
                    truths.append(False)
            new_mat[col] = truths

        self.corr_mat = new_mat


    def calc_max_ari(self, labels,
            n_neighbors=[2, 5, 10, 20, 50], def_res=None):
        """
        calculate the max ARI given organelle ground truths, and
        UMAP-embedding + leiden clustering.
        cycles through a given list of n_neighbors for the UMAP algorithm.

        Input labels is a dataframe and assumes merge key is in gene names
        """

        enrichments = self.enrichment_table.copy()
        enrichments.rename({'Gene names': 'gene_names'}, axis=1, inplace=True)

        metadata = enrichments['metadata'].copy()

        # merge organelle labels by finding all matching proteins in proteingroups
        organelles = []
        for _, row in metadata.iterrows():
            # combinations of gene names are attached by semicolons
            genes = str(row.gene_names).split(';')
            label_match = labels[labels['protein'].isin(genes)]
            if label_match.shape[0] > 0:
                # take (usually only) matching organelle
                organelles.append(label_match['organelle'].iloc[0])
            else:
                organelles.append('none')

        metadata['organelles'] = organelles

        enrichments.reset_index(drop=True, inplace=True)
        metadata.reset_index(drop=True, inplace=True)

        adata = ad.AnnData(enrichments['sample'], obs=enrichments['metadata'])
        aris = []
        resolutions = []
        for n in n_neighbors:
            output = clustering_workflows.calculate_max_ari(
                adata, n, 'organelles', res=0.4, n_random_states=2, def_res=def_res)
            aris.append(output[0])
            resolutions.append(output[1])
        max = np.max(aris)
        max_idx = aris.index(max)
        res = resolutions[max_idx]
        n_neighbor = n_neighbors[max_idx]

        return max, res, n_neighbor


    def grouped_reference_testing(self, labels, condition='-infected',
            merge_col='Gene names', label_col='organelle'):
        """
        Using the enrichment table, test the organellar difference between
        a contrast and control group, return a pval/enrichment table.
        """

        enrichment = self.enrichment_table.copy()
        labels = labels[[merge_col, label_col]].copy()

        # parse samples that are experiment controls
        samples = list(enrichment['sample'])
        control_samples = [x for x in samples if condition not in x]
        condition_samples = [x for x in samples if condition in x]
        control_samples.sort()

        # control-experiment pair dictionary
        sample_dict = {}
        for sample in control_samples:
            condition_pair = [x for x in condition_samples if sample in x][0]
            sample_dict[sample] = condition_pair

        # merge enrichment table with labels
        enrichment = enrichment.droplevel(0, axis=1)
        merged = enrichment.merge(labels, on=merge_col, how='left')

        # get unique labels except nans
        orgs = merged[label_col].unique[1:]
        pvals = pd.DataFrame()
        pvals['organelle'] = orgs

        # find t-test significance by each organelle
        for sample in sample_dict.keys():
            cond_sample = sample_dict[sample]

            control = merged[[sample, label_col]].copy()
            conditioned = merged[[cond_sample, label_col]].copy()

            pvs = []
            for org in orgs:
                control_orgs = control[control[label_col] == org]
                condition_orgs = conditioned[conditioned[label_col] == org]
                pval = scipy.stats.ttest_ind(
                    control_orgs[sample], condition_orgs[cond_sample])[1]
                pval = np.round(-1 * np.log10(pval), 2)
                pvs.append(pval)

            pvals[sample] = pvs

        self.ref_pvals_table = pvals























class RForestScoring():
    """
    This is a class for scoring of a dataset based on default RandomForest
    prediction algorithm, practically testing how much information is in the
    dataset for prediction power.

    """
    def __init__(self, dataset, reference_col, quant_cols):
        self.table = dataset
        self.ref_col = reference_col
        self.quant_cols = quant_cols


    def multiclass_balance_scale_split_predict(self, max_sample=150):
        # initiate variables
        table = self.table.copy()
        ref_col = self.ref_col
        quant_cols = self.quant_cols.copy()

        all_cols = quant_cols + [ref_col]
        table = table[all_cols]

        # find value counts of a label
        table.dropna(inplace=True)
        num_labels = table[ref_col].nunique()
        refs = table[ref_col].unique()
        refs.sort()

        # limit samples to max number alloted, this is semi-balancing
        samples = []
        for ref in refs:
            sample = table[table[ref_col] == ref]
            if sample.shape[0] <= max_sample:
                samples.append(sample)
            else:
                random_sample = sample.sample(max_sample)
                samples.append(random_sample)

        # concat all the reference samples
        balanced = pd.concat(samples).reset_index(drop=True)

        labels = pd.factorize(balanced[ref_col])
        definitions = labels[1]

        balanced['label'] = labels[0]

        # scale the quant columns with a StandardScaler
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


    def one_balance_scale_split_predict(self, ref):
        # initiate variables
        table = self.table.copy()
        ref_col = self.ref_col
        quant_cols = self.quant_cols.copy()

        table.dropna(inplace=True)
        refs = table[ref_col].unique()
        refs.sort()

        # count sample size for the reference being tested
        sample = table[table[ref_col] == ref]
        sample_size = sample.shape[0]

        # balance reference with all other samples
        others = table[table[ref_col] != ref]
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


    def repeat_collect_tests(self, one_vs_all=False, ref=None, max_sample=150, repeats=100):
        """
        Repeat either 1 v all or multiclass random forest tests and return
        confusion tables for precision and recall
        """
        tests = []
        predictions = []

        # repeat random forest prediction tests and save prediction results
        for i in np.arange(repeats):
            if one_vs_all:
                y_tested, y_predicted = self.one_balance_scale_split_predict(ref=ref)
            else:
                y_tested, y_predicted = self.multiclass_balance_scale_split_predict(
                    max_sample=max_sample)
            tests.append(y_tested)
            predictions.append(y_predicted)

        # concatenate all prediction results
        all_tests = np.concatenate(tests)
        all_preds = np.concatenate(predictions)

        # generate confusion matrix
        recall_table, precision_table = self.confusion_precision_recall(
            all_tests, all_preds)

        return recall_table, precision_table


    def confusion_precision_recall(tests, predictions, exp=''):
        """
        class method to create confusion chart from the output of repeat_collect_tests
        function
        """
        cross_recall = pd.crosstab(
            tests, predictions, rownames=['Actual Compartment'],
            colnames=['Predicted Compartment']).apply(
                lambda r: np.round(r/r.sum(), 3), axis=1)

        cross_precision = pd.crosstab(
            tests, predictions, rownames=['Actual Compartment'],
            colnames=['Predicted Compartment']).apply(
                lambda r: np.round(r/r.sum(), 3), axis=0)

        orgs = list(cross_recall)
        recall_table = {}
        precision_table = {}
        for org in orgs:
            recall_table[org] = cross_recall[org][org]
            precision_table[org] = cross_precision[org][org]

        recall_table = pd.DataFrame(recall_table, index=[exp])
        precision_table = pd.DataFrame(precision_table, index=[exp])

        return recall_table, precision_table
