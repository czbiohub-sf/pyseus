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

import primary_analysis as pa

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



class Validation():
    """
    Validation class takes as input standard_hits_table from
    AnalysisTables class (from primary_analysis.py) as input tables for various 
    post-processing or validation methods
    """

    def __init__(self, hit_table, corum, localization_table):
        """
        initiate class with a hit table (without interactions called) and other 
        necessary tables required for precision-recall analysis =
        """

        self.hit_table = hit_table
        self.corum = corum
        self.localization_table = localization_table
    
    def static_fdr(self, curvature, offset):
        """
        Call significant interactors from standard hits table with user input
        offset and curvature for thresholding
        """

        hits = self.hit_table.copy()
        hits['fdr'] = None
        hits['fdr'] = hits['fdr'].apply(lambda x: [curvature, offset])

        bait_pval = hits['pvals']
        enrichment = hits['enrichment']

        threshold = enrichment.apply(pa.calc_thresh, args=[curvature, offset])
        hits['interaction'] = np.where(
            (bait_pval > threshold) > True, False)

        self.fixed_interaction_table = hits

    def dynamic_fdr(self, perc=10, curvature=3, offset_seed=2.5):
        """
        compute dynamic FDR for each plate-bait group in all_hits table
        all_hits: DataFrame, output of hit_calling_validations.get_all_interactors
        perc: threshold for get_fdr5_seed specifying how much % of false positives
            are allowed before calling a FDR threshold
        RETURN:
        DataFrame, dataframe containing all the major and minor FDR thresholds
            for each bait-plate group
        DataFrame, all_hits df with fdr threshold columns added
        """
        hits = self.hit_table.copy()
        
        # group hits table by experiment & target and place them into a bin of lists
        selects = []
        grouped = hits.groupby(['experiment', 'target'])
        baits = []
        experiments = []
        for bait, group in grouped:
            experiments.append(bait[0])
            baits.append(bait[1])
            selects.append(group)
        
        # parallel processing for calculating FDR seed

        p = Pool()
        seeds = p.starmap(dfdr_find_thresh, zip(selects, baits,
            repeat(perc), repeat(curvature), repeat(offset_seed)))
        p.close()
        p.join()

        fdr_full = [[curvature, seed] for seed in seeds]
        
        fdr_df = pd.DataFrame()
        fdr_df['experiment'] = experiments
        fdr_df['target'] = baits
        fdr_df['fdr'] = fdr_full
        fdr_df.set_index('bait', inplace=True)    
        
        ############ Have to reconcile proper hits table format #######
        new_groups = []
        for bait, group in grouped:
            group = group.copy()
            experiment = bait[0]
            target = bait[1]
            fdr = fdr_df[(fdr_df['experiment']==experiment) &
                (fdr_df['target']==target)].fdr.item()           

            
            bait_pval = group['pvals']
            enrichment = group['enrichment']

            thresh = enrichment.apply(pa.calc_thresh,
                args=[fdr[0], fdr[1]])

            group['interaction'] = np.where(
                (bait_pval > thresh), True, False)

            group['fdr'] = [fdr] * group.shape[0]
            new_groups.append(group)

        interaction_table = pd.concat(new_groups)
        self.dynamic_fdr_table = fdr_df
        self.interaction_table = interaction_table


def dfdr_find_thresh(select, bait, perc=10, curvature=3, seed=2.5):
    """
    Find the proper p-val/enrichment threshold for a bait
    """
    
    # filter for negative hits
    neg_select = select[select['enrichment'] < 0]
    pos_select = select[select['enrichment'] > 0]

    # calcuate initial hit count by given curvature and seed
    hit = hit_count(neg_select, curvature, seed)
    
    # Find a threshold that lies just outside one hit detection
    if hit > 0:
        while hit > 0 and seed < 10:
            if seed > 4.2:
                seed += 0.2
            else:
                seed += 0.1
            hit = hit_count(neg_select, 3, seed)
    else:
        while hit == 0:
            seed -= 0.1
            hit = hit_count(neg_select, 3, seed)
        seed += 0.1
    
    # With the calculated seed, find the threshold that meets the 
    # requirement of less than 2 neg hits or less than designated % of positive hits
    neg_hit = hit_count(neg_select, 3, seed)
    pos_hit = hit_count(pos_select, 3, seed)

    pos_perc = 100 * neg_hit / pos_hit

    while (neg_hit < 2 or pos_perc < perc) and seed > 0.1:

        seed -= 0.1
        neg_hit = hit_count(neg_select, 3, seed)
        pos_hit = hit_count(pos_select, 3, seed)
        pos_perc = 100 * neg_hit / pos_hit

    if pos_perc > perc:
        seed += 0.1

    return round(seed, 2)    



def hit_count(bait_series, curvature, offset):
    """
    Count # of hits possible in a bait series with a given curvature and offset
    """
    bait_series = bait_series.copy()
    thresh = bait_series['enrichment'].apply(calc_thresh, args=[curvature, offset])
    hit = np.where(bait_series['pvals'] > thresh, True, False)

    return hit.sum()

def calc_thresh(enrich, curvature, offset):
    """simple function to get FCD thresh to recognize hits"""
    if abs(enrich) < offset:
        return np.inf
    else:
        return curvature / (abs(enrich) - offset)
