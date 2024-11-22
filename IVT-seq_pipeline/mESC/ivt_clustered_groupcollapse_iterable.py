#!/usr/bin/env python3
import sys
import re
import gzip
import pandas as pd
import numpy as np
from collections import Counter
import glob
from fuzzywuzzy import fuzz
import Levenshtein

threshold_distance = 6 #max Levenshtein distance between barcodes at any position to be collapsed into a cluster

# Function to check if two strings are within the threshold distance
def is_within_distance(str1, str2):
    return fuzz.ratio(str1, str2) >= 100 - threshold_distance

def collapse_group(subset_df):
    
    group_collapsed = {}
    
    # Calculate count of unique 'read' values per 'umi' and add as new column
    group_collapsed['chrom'] = subset_df['chrom'].iloc[0]
    group_collapsed['pos'] = subset_df['pos'].iloc[0]
    
    mode_value_1 = subset_df['barcode1'].mode().values[0]
    mode_frequency = subset_df['barcode1'].value_counts(normalize=True)[mode_value_1]
    group_collapsed['barcode1'] = mode_value_1
    group_collapsed['barcode1_freq'] = mode_frequency
    
    mode_value_2 = subset_df['barcode2'].mode().values[0]
    mode_frequency = subset_df['barcode2'].value_counts(normalize=True)[mode_value_2] 
    group_collapsed['barcode2'] = mode_value_2
    group_collapsed['barcode2_freq'] = mode_frequency
        
    mode_value_3 = subset_df['strand_x'].mode().values[0]
    mode_frequency = subset_df['strand_x'].value_counts(normalize=True)[mode_value_3]
    group_collapsed['strand_x'] = mode_value_3
    group_collapsed['strand_x_freq'] = mode_frequency
        
    mode_value_4 = subset_df['strand_y'].mode().values[0]
    mode_frequency = subset_df['strand_y'].value_counts(normalize=True)[mode_value_4]
    group_collapsed['strand_y'] = mode_value_4
    group_collapsed['strand_y_freq'] = mode_frequency
    
    mode_value_5 = subset_df['combo_bc'].mode().values[0]
    mode_frequency = subset_df['combo_bc'].value_counts(normalize=True)[mode_value_5] 
    group_collapsed['combo_bc'] = mode_value_5
    group_collapsed['combo_bc_freq'] = mode_frequency
    
    mode_value_6 = subset_df['recombinase_site'].mode().values[0]
    mode_frequency = subset_df['recombinase_site'].value_counts(normalize=True)[mode_value_6] 
    group_collapsed['recombinase_site'] = mode_value_6
    group_collapsed['recombinase_site_freq'] = mode_frequency
    
    #checking that combination of barcodes equals most common barcode1 and barcode2
    combo_bc1 = mode_value_5[:20]
    combo_bc2 = mode_value_5[20:]
    
    if (combo_bc1 == mode_value_1) & (combo_bc2 == mode_value_2):
        group_collapsed['match'] = 'yes'
    else:
        group_collapsed['match'] = 'no'
    
    #counting all reads and umis from the cluster
    group_collapsed['read_count'] = subset_df['read'].nunique()
    group_collapsed['umi_count_read1'] = subset_df['umi_y'].nunique() #number of unique RNA molecules coming from PCR1. No correction
    group_collapsed['umi_count_read2'] = subset_df['umi_x'].nunique() #number of unique RNA molecules coming from start of RT. No correction
    group_collapsed['unique_len_count'] = subset_df['len_aligned'].nunique() #reads of unique lengths - similar to umi_count_read2
    
    if len(subset_df['allele'].unique()) == 1:
        # If only one unique value in 'allele' column, collapse the group to a single row
        group_collapsed['true_allele'] = subset_df['allele'].values[0]
        
        return pd.DataFrame(group_collapsed,index=[0])
    
    else:
        valid_values = subset_df[(subset_df['allele'] != 'noVariant') & (subset_df['allele'] != 'neither')]
        if len(valid_values) > 0:
            # Get the value with the highest frequency
            true_allele = valid_values['allele'].value_counts().idxmax()
            group_collapsed['true_allele'] = true_allele
            
            return pd.DataFrame(group_collapsed,index=[0])
        
        else:
            group_collapsed['true_allele'] = 'neither'
            return pd.DataFrame(group_collapsed,index=[0])

        
def cluster_groups(group):
    
    # Initialize clusters
    clusters = []

    # Iterate through rows in the DataFrame
    for index, row in group.iterrows():
        found = False
        for cluster in clusters:
            # Check if the current string is within the threshold distance of any string in the cluster
            if any(is_within_distance(row['combo_bc'], cluster_string) for cluster_string in cluster):
                cluster.append(row['combo_bc'])
                found = True
                break
        if not found:
            # Create a new cluster for the current string
            clusters.append([row['combo_bc']])
    
    clustered_df = pd.DataFrame() #initialize dataframe for all clustrs at position

    for i, cluster in enumerate(clusters):
        
        subset_df = group[group['combo_bc'].isin(cluster)]
        frames = [clustered_df,collapse_group(subset_df)]
        clustered_df = pd.concat(frames)
      
    return clustered_df

if __name__ == "__main__":
    variants = sys.argv[1]
    bc_df_name = variants.split('_ITR')[0].split('mapped/')[1]

    if variants.endswith('.gz'):
        # Read the gzip-compressed table using pandas
        mapped_df = pd.read_table(f'{variants}', compression='gzip')
        bc_df = pd.read_table((glob.glob(f"{bc_df_name}_*_bc.tsv"))[0], compression='gzip')

    else:
        mapped_df = pd.read_table(f'{variants}')    
        bc_df = pd.read_table((glob.glob(f"{bc_df_name}_*_bc.tsv"))[0])

    mapped_df['len_aligned'] = mapped_df['end'] - mapped_df['start']
    
    #removing all rows that contain neither - very small proportion 
    mapped_df = mapped_df[mapped_df.allele != 'neither']

    #removing @ symbol from bc df
    bc_df['identifier'] = bc_df['identifier'].str[1:]

    #merging mapped df (read 2) genomic locations and bc df (containing reads)
    merged_df = pd.merge(mapped_df,bc_df,left_on='read',right_on='identifier')

    #making sure umi is correct sequence.
    merged_df['umi_x'] = merged_df['seq_x'].str[0:8]

    merged_df['combo_bc']=merged_df['barcode1']+merged_df['barcode2']

    group_collapsed = pd.concat([cluster_groups(group) for _,group in merged_df.groupby(['chrom','pos'])]).sort_values(by=['chrom', 'pos', 'read_count'], ascending=[True, True, False]).reset_index()
    group_collapsed.to_csv(f'{bc_df_name}_group_collapsed_clustered.tsv',sep = '\t')
