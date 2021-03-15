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
import pval_calculation as pval
import hack_pvals as hp
import pandas as pd
import numpy as np
import pyseus as pys
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



def all_hits_get_fdrs(pvals, perc=10):
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

    pvals = pvals.copy()
    # baits = list(pvals['target'].unique())
    selects = []
    grouped = pvals.groupby('target')
    baits = []
    for bait, group in grouped:
        baits.append(bait)
        selects.append(group)

    print("starting FDR1 calc")

    p = Pool()
    fdr1_seeds = p.starmap(get_fdr1_seed, zip(selects, baits))
    p.close()
    p.join()

    print("FDR1 seed finished!")

    fdr1_full = [[3, seed] for seed in fdr1_seeds]

    fdr5_args = zip(selects, fdr1_seeds, repeat(perc))

    p = Pool()
    fdr5_seeds = p.starmap(get_fdr5_seed, fdr5_args)
    p.close()
    p.join()

    print("FDR5 calculated!")

    fdr5_full = [[3, seed] for seed in fdr5_seeds]

    fdr_df = pd.DataFrame()
    fdr_df['bait'] = baits
    fdr_df['fdr1'] = fdr1_full
    fdr_df['fdr5'] = fdr5_full
    fdr_df.set_index('bait', inplace=True)

    # Find hits for FDR1 and FDR2
    for i, bait in enumerate(baits):
        fdr1 = fdr_df.loc[bait]['fdr1']
        fdr5 = fdr_df.loc[bait]['fdr5']

        select = selects[i]
        bait_pval = select['pvals']
        enrichment = select['enrichment']

        # 1% thresh
        first_thresh = enrichment.apply(hp.calc_thresh,
            args=[fdr1[0], fdr1[1]])

        # 5% thresh
        second_thresh = enrichment.apply(hp.calc_thresh,
            args=[fdr5[0], fdr5[1]])

        select['hits'] = np.where(
            (bait_pval > first_thresh), True, False)

        select['minor_hits'] = np.where(
            ((bait_pval < first_thresh) & (bait_pval > second_thresh)), True, False)

        select['fdr1'] = [fdr1] * select.shape[0]
        select['fdr5'] = [fdr5] * select.shape[0]
        selects[i] = select

    pvals = pd.concat(selects)

    return fdr_df, pvals


def all_interactors_set_hits(pvals, fdr1, fdr5):
    pvals = pvals.copy()

    pvals['fdr1'] = pvals['fdr1'].apply(lambda x: fdr1)
    pvals['fdr5'] = pvals['fdr5'].apply(lambda x: fdr5)


    bait_pval = pvals['pvals']
    enrichment = pvals['enrichment']

    # 1% thresh
    first_thresh = enrichment.apply(hp.calc_thresh,
        args=[fdr1[0], fdr1[1]])

    # 5% thresh
    second_thresh = enrichment.apply(hp.calc_thresh,
        args=[fdr5[0], fdr5[1]])

    pvals['hits'] = np.where(
        (bait_pval > first_thresh), True, False)

    pvals['minor_hits'] = np.where(
        ((bait_pval < first_thresh) & (bait_pval > second_thresh)), True, False)

    return pvals


def bait_set_hits(pvals, bait, plate, fdr1, fdr5):
    pvals = pvals.copy()

    section = pvals[(pvals['plate'] == plate) & (pvals['target'] == bait)]
    rest = pvals[(pvals['plate'] != plate) | (pvals['target'] != bait)]


    section['fdr1'] = section['fdr1'].apply(lambda x: fdr1)
    section['fdr5'] = section['fdr5'].apply(lambda x: fdr5)


    bait_pval = section['pvals']
    enrichment = section['enrichment']

    # 1% thresh
    first_thresh = enrichment.apply(hp.calc_thresh,
        args=[fdr1[0], fdr1[1]])

    # 5% thresh
    second_thresh = enrichment.apply(hp.calc_thresh,
        args=[fdr5[0], fdr5[1]])

    section['hits'] = np.where(
        (bait_pval > first_thresh), True, False)

    section['minor_hits'] = np.where(
        ((bait_pval < first_thresh) & (bait_pval > second_thresh)), True, False)

    pvals = pd.concat([rest, section])

    return pvals


def get_fdrs(pvals):
    """get dynamic FDR values and recompute new hits"""

    pvals = pvals.copy()

    baits = pvals.columns.get_level_values('baits')
    baits = list(set(baits))
    baits = [bait for bait in baits if bait != 'gene_names']
    selects = []
    for bait in baits:
        selects.append(pvals[bait])
    p = Pool()
    fdr1_seeds = p.starmap(get_fdr1_seed, zip(selects, baits))
    p.close()
    p.join()

    print("FDR1 seed finished!")

    fdr1_full = [[3, seed] for seed in fdr1_seeds]

    fdr5_args = zip(selects, fdr1_seeds)

    p = Pool()
    fdr5_seeds = p.starmap(get_fdr5_seed, fdr5_args)
    p.close()
    p.join()

    print("FDR5 calculated!")

    fdr5_full = [[3, seed] for seed in fdr5_seeds]

    fdr_df = pd.DataFrame()
    fdr_df['bait'] = baits
    fdr_df['fdr1'] = fdr1_full
    fdr_df['fdr5'] = fdr5_full
    fdr_df.set_index('bait', inplace=True)

    # Find hits for FDR1 and FDR2
    for bait in baits:
        fdr1 = fdr_df.loc[bait]['fdr1']
        fdr5 = fdr_df.loc[bait]['fdr5']

        bait_pval = pvals[bait]['pvals']
        enrichment = pvals[bait]['enrichment']

        # 1% thresh
        first_thresh = enrichment.apply(hp.calc_thresh,
            args=[fdr1[0], fdr1[1]])

        # 5% thresh
        second_thresh = enrichment.apply(hp.calc_thresh,
            args=[fdr5[0], fdr5[1]])

        pvals[(bait, 'hits')] = np.where(
            (bait_pval > first_thresh), True, False)

        pvals[(bait, 'minor_hits')] = np.where(
            ((bait_pval < first_thresh) & (bait_pval > second_thresh)), True, False)

        pvals[(bait, 'fdr1')] = [fdr1] * pvals.shape[0]
        pvals[(bait, 'fdr5')] = [fdr5] * pvals.shape[0]


    pvals.sort_index(axis=1, inplace=True)

    return fdr_df, pvals


def get_fdr1_seed(select, bait):
    # print(bait)
    neg_select = select[select['enrichment'] < 0]
    # neg_select = neg_select[neg_select['pvals'] < 15]
    # neg_select = neg_select[neg_select['enrichment'] > -4.5]
    seed = 2.5

    hit = hit_count(neg_select, 3, seed)

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
    return round(seed, 2)


def get_fdr5_seed(select, fdr1_seed, perc=10):
    neg_select = select[select['enrichment'] < 0]
    # neg_select = neg_select[neg_select['pvals'] < 15]
    # neg_select = neg_select[neg_select['enrichment'] > -4.5]

    pos_select = select[select['enrichment'] > 0]


    seed = fdr1_seed

    neg_hit = hit_count(neg_select, 3, seed)
    pos_hit = hit_count(pos_select, 3, seed)

    pos_perc = 100 * neg_hit / pos_hit

    # while (neg_hit < 2 or pos_perc < 10) and seed > 0.1:
    while (neg_hit < 2 or pos_perc < perc) and seed > 0.1:

        seed -= 0.1
        neg_hit = hit_count(neg_select, 3, seed)
        pos_hit = hit_count(pos_select, 3, seed)
        pos_perc = 100 * neg_hit / pos_hit

    if pos_perc > perc:
        seed += 0.1

    return round(seed, 2)


def hit_count(bait_series, fdr1, fdr2):
    bait_series = bait_series.copy()
    thresh = bait_series['enrichment'].apply(hp.calc_thresh_negs, args=[fdr1, fdr2])
    hit = np.where(bait_series['pvals'] > thresh, True, False)

    return hit.sum()
