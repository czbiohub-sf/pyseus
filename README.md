# pyseus
### Pyseus: Perseus in Python


##### This repository contains source code and notebooks to process and analyze Mass-Spec output from MaxQuant software, mainly utilizing the proteingroups.txt output file. The protein groups file is converted into Pandas DataFrames stored as objects in appropriate python Classes, which contain transformation and analysis methods. There are currently 4 main classes, with brief description of each class. For complete description and required arguments, please refer to the source code:

### RawTables
RawTables class takes the initial proteingroups.txt file and processes it into an organized, filtered DataFrame ready for significance/enrichment calculations. The following are some of the main class methods:

* filter_table() : filter rows that do not meet the QC (contaminants, reverse seq, only identified by site), and also filter non-intensity columns that will not be used for further processing
* transform_intensities(): transform intensities values by any given function, default is log2 transformation
* group_replicates(): Using regular expression, group replicates of the same biological sample together through multi-level columns in a dataframe
* remove_invalid_rows(): Remove rows (protein groups) that do not have at least one group that has values in all triplicates
* bait_impute(): bait-imputation for sets of data without enough samples- imputes a value from a normal distribution of the left-tail of a bait’s capture distribution for the undetected preys using multi-processing.
* prey_impute(): For protein groups with less than threshold number of sample number, impute a value from a normal distribution of the prey’s capture distribution using multi-processing. Note- most protein groups do not need imputation with 12-plate MBR

### AnalysisTables
AnalysisTables class uses imputed table and exclusion matrix output from RawTables, and performs essential analysis including enrichment/significance testing to call interactions and stoichiometry calculations. 


##### Exclusion Matrix manipulation: AnalysisTables class can utilize exclusion matrix, which is a user-defined set of sample & negative control to be used for enrichment/significance testing. The following are helpful methods:
* restore_default_exclusion_matrix(): restore the exclusion matrix to default - no exclusion of any baits / negative controls
* load_exclusion_matrix(): load user-defined exclusion matrix
* print_baits(): print all baits in the analysis
* print_controls(): print all the selected controls for an input bait
* print_excluded_controls(): print all the excluded controls for an input bait
* select_wildtype_controls(): Using string operation, select only matching samples to use as controls and exclude all others. Default search string: '_WT'

##### enrichment / signifiance testing
* simple_pval_enrichment(): calculate enrichment and pvals for each bait using user-defined exclusion matrix
* two_step_bootstrap_pval_enrichment: The two-step bootstrap pval/enrichment calculations does not use an exclusion table of user defined controls. It automatically drops statistically significant outliers from the prey pool on the first round, and uses the distribution with dropped outliers to calculate bootstrapped null distribution. The null distribution of the preys are then used in the second round for pval and enrichment calculation.
* convert_to_standard_table(): the standard table no longer uses column organization for baits. It follows a more SQL-like form where bait information is provided in separate columns

### Validation
Validation class takes standard_hits_table from AnalysisTables class for various post-processing or validation methods
* static_fdr(): Call significant interactors from standard hits table with user input offset and curvature for thresholding
* dynamic_fd(): compute dynamic FDR for each plate-bait group in the standard table


## Included Notebooks
Notebooks are useful for streamlined analysis, often taking wrapper functions that encompass many methods together. For proper usage, the user must save directories that include the proteingroups.txt file in the data folder. 

#### change_sample_names.ipynb 
This notebook allows easy manipulation of column names before raw processing. Pyseus works with intensity columns that are organized in <b>experiment_sample_replicate# </b>. For example <b>Glycine-Low-pH_LAMP1_2</b> would be organized as '<b>Glycine-Low-pH</b>' (<i>experiment</i>), '<b>LAMP1</b>' (<i>sample</i>), and '<b>2</b>' (<i>replicate #</i>). Dashes can be used as added descriptors for each section, but underscores need to be reserved for the separation. ex: 'Glycine_Low_pH_LAMP1_2 will not work as a sample title. 
#### process_raw_file

