import urllib.parse
import urllib.request
import sys
import multiprocessing
import os
import re
import textwrap
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import Pool


def czb_initial_processing(file_name, root, analysis, intensity_type='LFQ intensity',
    bait_impute=True, distance=1.8, width=0.3, thresh=100):
    
    """
    wrapper script for all the pre-processing up to imputation using
    PyseusRawTables Class. Saves and returns the PyseusRawTables in the
    designated analysis directory
    """
    # make directory for analysis folder
    analysis_dir = root + analysis

    if not os.path.isdir(analysis_dir):
        os.mkdir(analysis_dir)
    
    # Run all the processing methods
    pyseus_tables = PyseusRawTables(file_name=file_name,
        analysis=analysis, intensity_type=intensity_type)
    pyseus_tables.filter_table()
    pyseus_tables.transform_intensities(func=np.log2)
    pyseus_tables.group_replicates(intensity_re=r'_\d+$', reg_exp=r'(.*_.*)_\d+$')
    pyseus_tables.remove_invalid_rows()
    if bait_impute:
        pyseus_tables.bait_impute(distance=distance, width=width)
    else:
        pyseus_tables.prey_impute(distance=distance, width=width, thresh=thresh)
    
    # Save the class file in the analysis folder
    


class PyseusRawTables:
    """
    The Raw Tables class contain DataFrame objects and functions that cover
    multiple pre-processing steps to create a final processed imputed table. 
    
    """
    
    # initiate raw table by importing from data directory
    def __init__(self, file_name, analysis, intensity_type):
        
        self.raw_table = pd.read_csv(file_name,
            sep='\t', index_col=0, header=0, low_memory=False)
        # analysis string
        self.analysis = analysis
        # Specificastion of which intensity (raw or LFQ) to use
        self.intensity_type = intensity_type
      
    def filter_table(self, verbose=True):
        """filter rows that do not meet the QC (contaminants, reverse seq, only identified by site)
        Also filter non-intensity columns that will not be used for further processing"""
        
        ms_table = self.raw_table.copy()

        pre_filter = ms_table.shape[0]

        # remove rows with potential contaminants
        ms_table = ms_table[ms_table['Potential contaminant'].isna()]

        # remove rows only identified by site
        ms_table = ms_table[ms_table['Only identified by site'].isna()]

        # remove rows that are reverse seq
        ms_table = ms_table[ms_table['Reverse'].isna()]

        filtered = pre_filter - ms_table.shape[0]
        if verbose:
            print("Filtered " + str(filtered) + ' of '
                  + str(pre_filter) + ' rows. Now '
                  + str(ms_table.shape[0]) + ' rows.')
        
        # select necessary columns
        all_cols = list(ms_table)
        info_cols = ['Protein IDs', 'Majority protein IDs', 'Protein names',
        'Gene names', 'Fasta headers']
        intensity_cols = select_intensity_cols(all_cols, self.intensity_type)

        ms_table = ms_table[info_cols + intensity_cols]


        self.filtered_table = ms_table
    
    def transform_intensities(self, func=np.log2):
        """transform intensity values in the dataframe to a given function"""
        
        try:
            filtered = self.filtered_table.copy()
        except AttributeError:
            print(
                "Raw table has not been filtered yet, use filter_table() method"\
                "before transforming intensities")
            return

        filtered = self.filtered_table.copy()
        intensity_cols = select_intensity_cols(list(self.filtered_table),
            intensity_type=self.intensity_type)

        # for each intensity column, transform the values
        for int_col in intensity_cols:
            # if transformation is log2, convert 0s to nans
            # (faster in one apply step than 2)
            if func == np.log2:
                filtered[int_col] = filtered[int_col].apply(lambda x: np.nan
                    if x == 0 else func(x))
            else:
                filtered[int_col] = filtered[int_col].apply(func)
                # Replace neg inf values is np.nan
                filtered[int_col] = filtered[int_col].apply(
                    lambda x: np.nan if np.isneginf(x) else x)
        
        self.transformed_table = filtered
    
    def group_replicates(self, intensity_re=r'_\d+$', reg_exp=r'(.*_.*)_\d+$'):
        """Group the replicates of intensities into replicate groups"""
        
        reg_exp = self.intensity_type + reg_exp
        try: 
            self.transformed_table
            transformed = self.transformed_table.copy()
        except AttributeError:
            print(
            'Intensity values have not been transformed yet from '\
            'filtered table,\nwe recommend using transform_intensities() '\
            'method before grouping replicates.\n')

            try: 
                self.filtered_table
                print("Using filtered_table to group replicates.")
                transformed = self.filtered_table.copy()
            except AttributeError:
                print('Please filter raw table first using filter_table()\
                    method.')
                return


       # get col names
        col_names = list(transformed)

        # using a dictionary, group col names into replicate groups
        group_dict = {}
        for col in col_names:
            # search REs of the replicate ID, and get the group names

            # search if the col is for intensity values
            intensity_search = re.search(intensity_re, col.lower(),
                flags=re.IGNORECASE)

            # if so, get the group name and add to the group dict
            # use groups from re.search to customize group names
            if intensity_search:
                group_search = re.search(reg_exp, col, flags=re.IGNORECASE)
                group_name = ''

                for re_group in group_search.groups():
                    group_name += re_group
                group_dict[col] = group_name

            # if not, group into 'Info'
            else:
                group_dict[col] = 'Info'


        # pd function to add the replicate group to the columns
        grouped = pd.concat(dict((*transformed.groupby(group_dict, 1),)), axis=1)

        grouped.columns = grouped.columns.rename("Baits", level=0)
        grouped.columns = grouped.columns.rename("Replicates", level=1) 

        self.grouped_table = grouped
    
    
    def remove_invalid_rows(self):
        """Remove rows that do not have at least one group that has values
        in all triplicates"""

        try:
            grouped = self.grouped_table.reset_index(drop=True).copy()
        except AttributeError:
            print("Replicates need to be grouped before this method."\
                "Please use group_replicates() to group replicates under same sample")
            return

        # reset index
        grouped = self.grouped_table.reset_index(drop=True).copy()
        unfiltered = self.grouped_table.shape[0]

        # Get a list of all groups in the df
        group_list = list(set([col[0] for col in list(grouped) if col[0] != 'Info']))

        # booleans for if there is a valid value
        filtered = grouped[group_list].apply(np.isnan)
        # loop through each group, and filter rows that have valid values
        for group in group_list:
            # filter all rows that qualify as all triplicates having values
            filtered = filtered[filtered[group].any(axis=1)]

        # a list containing all the rows to delete
        del_list = list(filtered.index)

        # create a new df, dropping rows with invalid data
        filtered_df = grouped.drop(del_list)
        filtered_df.reset_index(drop=True, inplace=True)
        filtered = filtered_df.shape[0]

        print("Removed invalid rows. " + str(filtered) + " from "
            + str(unfiltered) + " rows remaining.")

        self.preimpute_table = filtered_df
    
    def bait_impute(self, distance=1.8, width=0.3):
        """
        bait-imputation for sets of data without enough samples.
        This fx imputes a value from a normal distribution of the left-tail
        of a bait’s capture distribution for the undetected preys using
        multi-processing.
            distance: float, distance in standard deviation from the
            mean of the sample distribution upon which to impute. Default = 0
            width: float, width of the distribution to impute in standard deviations. Default = 0.3
        """
        
        try:
            imputed = self.preimpute_table.copy()
        except AttributeError:
            print("group_replicates() and remove_invalid_rows() need to be run"\
                "before imputation")
            return
        
        self.bait_impute_params = {'distance': distance, 'width': width}

        # Retrieve all col names that are not classified as Info
        bait_names = [col[0] for col in list(imputed) if col[0] != 'Info']
        baits = list(set(bait_names))
        bait_series = [imputed[bait].copy() for bait in baits]
        bait_params = zip(
            bait_series, repeat(distance), repeat(width))

        # Use multiprocessing pool to parallel impute
        p = Pool()
        impute_list = p.starmap(pool_impute, bait_params)
        p.close()
        p.join()

        for i, bait in enumerate(baits):
            imputed[bait] = impute_list[i]

        self.bait_imputed_table = imputed
    
    def prey_impute(self, distance=0, width=0.3, thresh=100):
        """
        default mode of imputation. For protein groups with less than threshold number
        of sample number, impute a value from a normal distribution of the prey’s capture
        distribution using multi-processing. Note- most protein groups do not need imputation
        with 12-plate MBR

            distance: float, distance in standard deviation from the mean of the
                sample distribution upon which to impute. Default = 0
            width: float, width of the distribution to impute in standard deviations.
                Default = 0.3
            threshold: int, max number of samples required for imputation
        """
        
        try:
            imputed = self.preimpute_table.copy()
        except AttributeError:
            print("group_replicates() and remove_invalid_rows() need to be run"\
                "before imputation")
            return
        
        imputed = self.preimpute_table.copy()
        imputed.drop(columns='Info', inplace=True)
        imputed = imputed.T
        self.prey_impute_params = {'distance': distance, 'width': width,
            'thresh': thresh}

        # Retrieve all col names that are not classified as Info
        baits = list(imputed)
        bait_series = [imputed[bait].copy() for bait in baits]
        bait_params = zip(
            bait_series, repeat(distance), repeat(width), repeat(thresh))

        # Use multiprocessing pool to parallel impute
        p = Pool()
        impute_list = p.starmap(pool_impute_prey, bait_params)
        p.close()
        p.join()

        for i, bait in enumerate(baits):
            imputed[bait] = impute_list[i]

        imputed = imputed.T

        info_cols = [x for x in list(self.preimpute_table) if x[0] == 'Info']
        for col in info_cols:
            imputed[col] = self.preimpute_table[col]      

        self.prey_imputed_table = imputed  
        
def select_intensity_cols(orig_cols, intensity_type):
    """from table column names, return a list of only intensity cols
    rtype: intensity_cols list """
    # new list of intensity cols
    intensity_cols = []

    # create a regular expression that can distinguish between
    # intensity and LFQ intensity
    re_intensity = '^' + intensity_type.lower()

    # for loop to include all the intensity col names
    intensity_type = intensity_type.lower()
    for col in orig_cols:
        col_l = col.lower()

        # check if col name has intensity str
        if re.search(re_intensity, col_l):
            intensity_cols.append(col)

    return intensity_cols


def pool_impute(bait_group, distance=1.8, width=0.3):
    """target for multiprocessing pool from multi_impute_nans"""
    all_vals = bait_group.stack()
    mean = all_vals.mean()
    stdev = all_vals.std()

    # get imputation distribution mean and stdev
    imp_mean = mean - distance * stdev
    imp_stdev = stdev * width

    # copy a df of the group to impute values
    bait_df = bait_group.copy()

    # loop through each column in the group
    for col in list(bait_df):
        bait_df[col] = bait_df[col].apply(random_imputation_val,
            args=(imp_mean, imp_stdev))
    return bait_df


def pool_impute_prey(bait_group, distance=0, width=0.3, thresh=100):
    """target for multiprocessing pool from multi_impute_nans"""

    if bait_group.count() > thresh:
        return bait_group


    mean = bait_group.mean()
    stdev = bait_group.std()

    # get imputation distribution mean and stdev
    imp_mean = mean - distance * stdev
    imp_stdev = stdev * width

    # copy a df of the group to impute values
    bait_df = bait_group.copy()


    bait_df = bait_df.apply(random_imputation_val,
            args=(imp_mean, imp_stdev))
    return bait_df


def random_imputation_val(x, mean, std):
    """from a normal distribution take a random sample if input is
    np.nan. For real values, round to 4th decimal digit.
    Floats with longer digits will be 'barcoded' by further digits

    rtype: float"""

    if np.isnan(x):
        return np.random.normal(mean, std, 1)[0]
    else:
        return np.round(x, 4)
