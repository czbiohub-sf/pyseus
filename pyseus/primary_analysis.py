import urllib.parse
import urllib.request
import sys
import pdb
import collections
import multiprocessing
import itertools
import scipy
import random
import re
import pandas as pd
import numpy as np
# from pyseus import basic_processing as pys
from multiprocessing import Queue
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from itertools import repeat
from multiprocessing import Pool
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.stats import percentileofscore
from sklearn.metrics.pairwise import cosine_similarity


class AnalysisTables:
    """
    Analysis Tables contains DataFrame objects, functions, and metadata that cover
    essential analysis including enrichment/significance testing to call interactions 
    and stoichiometry calculations
    """

    def __init__(
        self,
        root,
        analysis,
        imputed_table,
        exclusion_matrix):

        # initiate class that cover essential metadata and imputed table
        # from RawTables class.
         
        self.root = root
        self.analysis = analysis
        self.imputed_table = imputed_table
        self.exclusion_matrix = exclusion_matrix

    def restore_default_exclusion_matrix(self):
        """
        Restore exclusion matrix to default - No exclusion
        """
        
        exclusion = self.exclusion_matrix.copy()
        baits = list(exclusion)
        baits.remove('Baits')
        for bait in baits:
            exclusion[bait] = True

        self.exclusion_matrix = exclusion


    def load_exclusion_matrix(self, alt_text=''):
        """
        Load user-defined exclusion_matrix for simple analysis
        """
        exclusion_matrix = pd.read_csv(self.root + self.analysis +
            '/analysis_exclusion_matrix'+ alt_text + '.csv')
        
        self.exclusion_matrix = exclusion_matrix
    
    def print_baits(self):
        """
        Show dataframe of all possible baits
        """
        return self.exclusion_matrix['Baits']
    
    def print_controls(self, bait):
        """
        Print all the selected controls for an input bait
        """

        excluded = self.exclusion_matrix.copy()
        excluded = excluded[['Baits', bait]]
        excluded = excluded[excluded[bait] == True]
        
        if excluded.shape[0] > 0:
            return excluded[['Baits']]
        else:
            print("No control baits selected as control")
    
    def print_excluded_controls(self, bait):
        """
        Print all the excluded controls for an input bait
        """

        excluded = self.exclusion_matrix.copy()
        excluded = excluded[['Baits', bait]]
        excluded = excluded[excluded[bait] == False]
        
        if excluded.shape[0] > 0:
            return excluded[['Baits']]
        else:
            print("No excluded baits in control")
    
    def select_wildtype_controls(self, wt_re='_WT'):
        """
        Using string operation, select only wildtypes to use as controls
        and exclude all others. Since this is based on string,
        one can customize any label to use for control samples.

        Does not override default excluded controls.
        """

        exclusion = self.exclusion_matrix.copy()
        exclusion.set_index('Baits', inplace=True)
        exclusion = exclusion.T  
        baits = list(exclusion)

        for bait in baits:
            if wt_re in bait:
                continue
            else:
                exclusion[bait] = False
        
        exclusion = exclusion.T
        exclusion.reset_index(inplace=True)

        self.exclusion_matrix = exclusion


    def simple_pval_enrichment(self, std_enrich=True, mean=False):
        """
        Calculate enrichment and pvals for each bait, no automatic removal
        """
        imputed = self.imputed_table.copy()
        exclusion = self.exclusion_matrix.copy()

        # iterate through each cluster to generate neg con group
        bait_list = [col[0] for col in list(imputed) if col[0] != 'Info']
        bait_list = list(set(bait_list))
        total = len(bait_list)
        baitrange = list(np.arange(total))   
          
        multi_args = zip(bait_list, repeat(imputed), baitrange, repeat(total),
        repeat(exclusion), repeat(std_enrich), repeat(mean))

        p = Pool()
        print("P-val calculations..")
        outputs = p.starmap(simple_pval, multi_args)
        p.close()
        p.join()
        print("Finished!")  

        master_df = pd.concat(outputs, axis=1)

        # join gene names to the df
        gene_names = imputed[[('Info', 'Protein IDs'), ('Info', 'Gene names')]]
        gene_names.set_index(('Info', 'Protein IDs'), drop=True, inplace=True)
        gene_names.rename(columns={'Info': 'gene_names'}, inplace=True)
        gene_names.rename(columns={'Gene names': 'gene_names'}, inplace=True)


        master_df = pd.concat([master_df, gene_names], axis=1, join='inner')

        self.simple_pval_table = master_df
    

    def convert_to_standard_table(self, metrics=['pvals', 'enrichment'], interactors=False,
            simple_analysis=True):
        """
        the standard table no longer uses column organization for baits. 
        It follows a more SQL-like form where bait information is provided in 
        separate columns.
        """
        if simple_analysis:
            pvals = self.simple_pval_table.copy()
            protein_ids = pvals.index.to_list()

        pvals.set_index(('gene_names', 'gene_names'), inplace=True)
        pvals.index.name = 'gene_names'
        targets = pvals.columns.get_level_values('baits')
        targets = list(set(targets))

        all_hits = []
        # Get all hits and minor hits along with the metric data
        for target in targets:
            target_pvs = pvals[target]
            target_pvs['protein_ids'] = protein_ids

            # return target_pvs
            # just_hits bool will return all hits, else it will only return interactors
            
            if interactors:
                selection = ['interactors'] + metrics
                hits = target_pvs[target_pvs['hits'] | target_pvs['minor_hits']][selection]
                hits.reset_index(inplace=True)
            else:
                hits = target_pvs
                hits.reset_index(inplace=True)


            hits['target'] = target.upper()
            hits.rename(columns={'gene_names': 'prey'}, inplace=True)
            hits.reset_index(drop=True, inplace=True)
            all_hits.append(hits)

        all_hits = pd.concat(all_hits, axis=0)
        
        if interactors:
            self.standard_interactors_table = all_hits
        else:
            self.standard_hits_table = all_hits
     

def simple_pval(bait, df, num, total, exclusion, std_enrich=True, mean=False):
    """ A first round of pval calculations to remove any significant hits
    from negative controls """

    df = df.copy()
    excluded = exclusion.copy()

    # initiate other variables required for the fx
    gene_list = df[('Info', 'Protein IDs')].tolist()

    # construct a negative control
    temporary = df.copy()
    neg_control = df.copy()
    temporary.drop('Info', level='Baits', inplace=True, axis=1)
    neg_control.drop('Info', level='Baits', inplace=True, axis=1)
 
    # Get a list of excluded genes
    excluded = excluded[['Baits', bait]]
    excluded = excluded[excluded[bait] == False]

    if excluded.shape[0] > 0:
        exclude_list = excluded['Baits'].to_list()
      
        # Convert all values in same groups as np.nans
        for gene in exclude_list:
            neg_control[gene] = neg_control[gene].where(
                neg_control[gene] > 100, np.nan)

    # calculate the p-values

    # combine values of replicates into one list
    bait_series = temporary[bait].values.tolist()

    # add an index value to the list for locating neg_control indices
    for i in np.arange(len(bait_series)):
        bait_series[i].append(i)

    # perform the p value calculations, with bagging replacement for np.nans
    pval_series = pd.Series(bait_series, index=gene_list, name='pvals')

    pval_series = pval_series.apply(get_pvals, args=[neg_control.T, std_enrich, mean])

    pvals, enrichment = pval_series.apply(lambda x: x[0]), pval_series.apply(lambda x: x[1])
    pvals.name = 'pvals'
    enrichment.name = 'enrichment'

    # Find positive hits from enrichment and pval calculations
    pe_df = pd.concat([pvals, enrichment], axis=1)

    output = pd.concat([pe_df[['enrichment', 'pvals']]], keys=[bait],
        names=['baits', 'values'], axis=1)

    if num % 20 == 0:
        print(str(num) + ' / ' + str(total) + ' baits processed')

    return output
    


def get_pvals(x, control_df, std_enrich, mean=False):
    """This is an auxillary function to calculate p values
    that is used in enrichment_pval_dfs function

    rtype: pval float"""

    # get the index to access the right set of control intensities
    row = x[-1]
    neg_con = control_df[row].values.tolist()
    pval = scipy.stats.ttest_ind(x[:-1], control_df[row].values.tolist(),
    nan_policy='omit')[1]

    # negative log of the pvals
    pval = -1 * np.log10(pval)


    # calculate enrichment
    if std_enrich:
        std = np.nanstd(neg_con)
        if mean:
            enrichment = (np.nanmean(x[:-1]) - np.nanmean(neg_con)) / std
        else:
            enrichment = (np.nanmedian(x[:-1]) - np.nanmedian(neg_con)) / std


    else:
        if mean:
            enrichment = (np.nanmean(x[:-1]) - np.nanmean(neg_con))
        else:
            enrichment = (np.nanmedian(x[:-1]) - np.nanmedian(neg_con))

    return [pval, enrichment]
