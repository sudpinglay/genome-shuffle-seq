import gzip
import pandas as pd
import numpy as np


with open('samples.txt','r') as input_f:
	lines = input_f.readlines()

for i in lines:
	name = str(i.strip())
																												
	output_file = name+'_grouped.tsv'
	output_file_2 = name+'_grouped_merged.tsv'
	output_file_3 = name+'_grouped_merged_collapsed.tsv'

	df = pd.read_csv(name+'_extracted_bcs.tsv',sep='\t')
	df_group = df.groupby(['barcode1','barcode2','read1_cs','read2_cs','umi']).size().reset_index(name='readcount')

	#grouping based on barcode1 and barcode2, summing up all reads across all UMIs for each BC combo
	df_group_bc = df_group.groupby(['barcode1', 'barcode2','read1_cs','read2_cs'])['readcount'].sum().reset_index()
	df_group_bc.to_csv(output_file,index=False)

	#count number of unique umis per BC combo
	umi_group = df_group.groupby(['barcode1', 'barcode2','read1_cs','read2_cs'])['umi'].nunique().reset_index(name='umi_count')

	#merging both dfs to yield a df where each barcode combo is matched with total number of reads across all UMIs and total number of UMIs
	merged_df = pd.merge(df_group_bc, umi_group, on=['barcode1','barcode2','read1_cs','read2_cs'], how='inner')
	merged_df.to_csv(output_file_2,index=False)

	#merging merged_df to itself and then collapsing those rows where barcode1 in one row equals barcode 2 in another row and vice versa
	grouped_merged_df = merged_df.merge(merged_df, left_on='barcode1', right_on='barcode2', how='outer')
	grouped_merged_collapsed_df = grouped_merged_df[(grouped_merged_df['barcode1_x'] == grouped_merged_df['barcode2_y']) & (grouped_merged_df['barcode2_x'] == grouped_merged_df['barcode1_y'])]
	grouped_merged_collapsed_df.to_csv(output_file_3,index=False)