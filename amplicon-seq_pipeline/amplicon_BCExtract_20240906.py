#!/usr/bin/env python3
import sys
import re
import gzip
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import glob


def extract_sequences(read1_file, read2_file, output_file):
    
    with open(output_file, 'w') as file:
        file.write("umi\treal_barcode1\tread1_cs\tread1_recomb_site\treal_barcode2\tread2_cs\tread2_recomb_site\ttrue_recomb_site\n")
        
        i = 0
        num_reads = 0
        same_cs_count = 0
        attL_count = 0
        attR_count = 0

        with gzip.open(read1_file, 'rt') as file1, gzip.open(read2_file, 'rt') as file2:
            while True:
                # Read 4 lines at a time from each file
                lines1 = [file1.readline().strip() for _ in range(4)]
                lines2 = [file2.readline().strip() for _ in range(4)]

                # If we've reached the end of either file, break the loop
                if not lines1[0] or not lines2[0]:
                    break

                num_reads += 1

                barcode1 = None
                capture_seq_1 = None
                recombinase_site_1 = None
                barcode2 = None
                capture_seq_2 = None
                recombinase_site_2 = None
                cell_bc = None
                umi = None
                true_recomb_site = None
                real_barcode1 = None  
                real_barcode2 = None

                identifier_1 = lines1[0].split(' ')[0]
                sequence_1 = lines1[1].strip()
                identifier_2 = lines2[0].split(' ')[0]
                sequence_2 = lines2[1].strip()

                if identifier_1 == identifier_2:

                    ### Processing read1 ###
                    
                    if re.findall(r'TGAGC(.{20})GTGG', sequence_1[20:75]):
                        umi = sequence_1[:10]
                        capture_seq_1 = 'CS2'
                        barcode1 = re.findall(r'TGAGC(.{20})GTGG', sequence_1[20:75])[0]
                        recombinase_site_1 = 'attP'

                    elif re.findall(r'TGAGC(.{20})GGCC', sequence_1[20:75]):
                        umi = sequence_1[:10]
                        capture_seq_1 = 'CS2'
                        barcode1 = re.findall(r'TGAGC(.{20})GGCC', sequence_1[20:75])[0]
                        recombinase_site_1 = 'attB'

                    elif re.findall(r'AAAGC(.{20})TGGG', sequence_1[20:75]):
                        umi = sequence_1[:10]
                        capture_seq_1 = 'CS1'
                        barcode1 = re.findall(r'AAAGC(.{20})TGGG', sequence_1[20:75])[0]
                        recombinase_site_1 = 'attP'

                    elif re.findall(r'AAAGC(.{20})CCGG', sequence_1[20:75]):
                        umi = sequence_1[:10]
                        capture_seq_1 = 'CS1'
                        barcode1 = re.findall(r'AAAGC(.{20})CCGG', sequence_1[20:75])[0]
                        recombinase_site_1 = 'attB'

                    elif re.findall(r'TGAGC(.{20})ATAAC',sequence_1[20:75]):
                        umi = sequence_1[:10]
                        capture_seq_1 = 'CS2'
                        barcode1 = re.findall(r'TGAGC(.{20})ATAA',sequence_1[20:75])[0]
                        recombinase_site_1 = 'loxP'

                    elif re.findall(r'AAAGC(.{20})ATAAC',sequence_1[20:75]):
                        umi = sequence_1[:10]
                        capture_seq_1 = 'CS1'
                        barcode1 = re.findall(r'AAAGC(.{20})ATAA',sequence_1[20:75])[0]
                        recombinase_site_1 = 'loxP'

                    ### Processing read2 ###
                    if re.findall(r'AAAGC(.{20})CCGG', sequence_2[12:53]):
                        capture_seq_2 = 'CS1'
                        barcode2 = re.findall(r'AAAGC(.{20})CCGG', sequence_2[12:53])[0]
                        recombinase_site_2 = 'attB'

                    elif re.findall(r'AAAGC(.{20})TGGG', sequence_2[12:53]):
                        capture_seq_2 = 'CS1'
                        barcode2 = re.findall(r'AAAGC(.{20})TGGG', sequence_2[12:53])[0]
                        recombinase_site_2 = 'attP'

                    elif re.findall(r'TGAGC(.{20})GTGG', sequence_2[12:53]):
                        capture_seq_2 = 'CS2'
                        barcode2 = re.findall(r'TGAGC(.{20})GTGG', sequence_2[12:53])[0]
                        recombinase_site_2 = 'attP'

                    elif re.findall(r'TGAGC(.{20})GGCC', sequence_2[12:53]):
                        capture_seq_2 = 'CS2'
                        barcode2 = re.findall(r'TGAGC(.{20})GGCC', sequence_2[12:53])[0]
                        recombinase_site_2 = 'attB'

                    elif re.findall(r'TGAGC(.{20})ATAAC',sequence_2[12:53]):
                        capture_seq_2 = 'CS2'
                        barcode2 = re.findall(r'TGAGC(.{20})ATAA',sequence_2[12:53])[0]
                        recombinase_site_2 = 'loxP'

                    elif re.findall(r'AAAGC(.{20})ATAAC',sequence_2[12:53]):
                        capture_seq_2 = 'CS1'
                        barcode2 = re.findall(r'AAAGC(.{20})ATAA',sequence_2[12:53])[0]
                        recombinase_site_2 = 'loxP'

                    
                    if capture_seq_1 and capture_seq_2:
                        if capture_seq_1 == capture_seq_2:
                            same_cs_count += 1

                    #### Assigning real_barcode sequences #####
                    if capture_seq_1 == 'CS2' and capture_seq_2 == 'CS1':
                        real_barcode1 = barcode1
                        real_barcode2 = barcode2

                        #### Looking for attL or attR sequences ####
                        if recombinase_site_1 == recombinase_site_2:
                            true_recomb_site = recombinase_site_1
                        else:
                            if recombinase_site_1 == 'attB' and recombinase_site_2 == 'attP':
                                true_recomb_site = 'attL'
                                attL_count+=1

                            elif recombinase_site_1 == 'attP' and recombinase_site_2 == 'attB':
                                true_recomb_site = 'attR'
                                attR_count+=1

                    elif capture_seq_1 == 'CS1' and capture_seq_2 == 'CS2':
                        real_barcode1 = barcode2
                        real_barcode2 = barcode1

                        #### Looking for attL or attR sequences ####
                        if recombinase_site_1 == recombinase_site_2:
                            true_recomb_site = recombinase_site_1
                        else:
                            if recombinase_site_1 == 'attB' and recombinase_site_2 == 'attP':
                                true_recomb_site = 'attR'
                                attR_count+=1

                            elif recombinase_site_1 == 'attP' and recombinase_site_2 == 'attB':
                                true_recomb_site = 'attL'
                                attL_count+=1

                    #### Counting CS1-CS1 or CS2-CS2 reads if any #####
                    
                    # Writing to output file
                    if real_barcode1 and real_barcode2:
                        
                        file.write(f"{umi}\t{real_barcode1}\t{capture_seq_1}\t{recombinase_site_1}\t{real_barcode2}\t{capture_seq_2}\t{recombinase_site_2}\t{true_recomb_site}\n")
                        i += 1

        fraction_used = i / num_reads if num_reads > 0 else 0
        print(f'number of extracted BCs in {sample}:{i}, fraction reads useful: {fraction_used}', flush=True)
        print(f'number of same CS in {sample}:{same_cs_count}', flush=True)
        print(f'number of attL in {sample}:{attL_count}', flush=True)
        print(f'number of attR in {sample}:{attR_count}', flush=True)
        print('\n\n', flush=True)


def group_collapse(extracted_bcs):
    df = pd.read_table(extracted_bcs)

    #keeping only those real_barcode1, real_barcode2 combos that were detected from both strands
    both_cs = df.groupby(['real_barcode1','real_barcode2']).filter(lambda x:x['read1_cs'].nunique()>1)

    df_group = both_cs.groupby(['real_barcode1','real_barcode2','true_recomb_site']).size().reset_index(name='read_count')
    df_umi = both_cs.groupby(['real_barcode1','real_barcode2','true_recomb_site'])['umi'].nunique().reset_index(name='umi_count')
    df_collapsed = pd.merge(df_group,df_umi,on=['real_barcode1','real_barcode2','true_recomb_site'])
    df_collapsed.to_csv(f'{sample}group_collapsed.tsv',sep='\t',index=False)


if __name__ == "__main__":
    sample = sys.argv[1]
    path = sys.argv[2]
    output_folder = sys.argv[3]
    
    read1_file = glob.glob(path + sample + '*R1_*.fastq.gz')[0]
    read2_file = glob.glob(path + sample + '*R2_*.fastq.gz')[0]
                                                                                                                        
    output_file = sample+'extracted_bcs.tsv'

    extract_sequences(read1_file, read2_file, output_file)
    group_collapse(output_file)


