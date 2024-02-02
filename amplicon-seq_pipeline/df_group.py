import gzip
import pandas as pd
import numpy as np


with open('samples.txt','r') as input_f:
	lines = input_f.readlines()

for i in lines:
	i = str(i.strip())
	name =  i.split('500_')[1].split('_S')[0]
	#read1_file = i+'_R1_001.fastq.gz'
	#read2_file = i+'_R2_001.fastq.gz'																															
	output_file = i+'_grouped.tsv'
	output_file_2 = i+'_grouped_merged.tsv'

	x = pd.read_csv(i+'.tsv',sep='\t')
	df = x.groupby(['barcode1','barcode2','umi']).size().reset_index(name='readcount')

	#grouping based on barcode1 and barcode2, summing up all reads across all UMIs for each BC combo
	group = df.groupby(['barcode1', 'barcode2'])['readcount'].sum().reset_index()
	df.to_csv(output_file,index=False)

	#count number of unique umis per BC combo
	umi_group = df.groupby(['barcode1', 'barcode2'])['umi'].nunique().reset_index(name='umi_count')

	#merging both dfs to yield a df where each barcode combo is matched with total number of reads across all UMIs and total number of UMIs
	merged_df = pd.merge(group, umi_group, on=['barcode1','barcode2'], how='inner')
	merged_df.to_csv(output_file_2,index=False)
	
