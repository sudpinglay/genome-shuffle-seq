import gzip
import pandas as pd
import numpy as np

def extract_sequences(read1_file, read2_file, output_file):
    read1_dict = {}
    read1_umi_dict = {}
    read1_cs_dict = {}
    
    read2_dict = {}
    read2_cs_dict = {}

    cs2_cs1_count = 0
    cs1_cs2_count = 0
    cs2_cs2_count = 0
    cs1_cs1_count = 0

    # Read and store sequences from Read 1 based on conditions
    with gzip.open(read1_file, 'rt') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            identifier = lines[i].split(' ')[0]
            sequence = lines[i + 1].strip()

            if (sequence[28:32] == 'GAGC' and sequence[52:56] == 'ATAA' and sequence[10:14] == 'CCTT'):
                capture_seq = 'CS2'
                barcode1 = sequence[32:52]
                umi = sequence[:10]
                
                read1_umi_dict[identifier] = umi
                read1_dict[identifier] = barcode1
                read1_cs_dict[identifier] = capture_seq

            elif (sequence[28:32] == 'AAGC' and sequence[52:56] == 'ATAA' and sequence[10:14] == 'TTGC'):
                capture_seq = 'CS1'
                barcode1 = sequence[32:52]
                umi = sequence[:10]
                
                read1_umi_dict[identifier] = umi
                read1_dict[identifier] = barcode1 
                read1_cs_dict[identifier] = capture_seq
                
    
    # Read and store sequences from Read 2 based on conditions
    with gzip.open(read2_file, 'rt') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            identifier = lines[i].split(' ')[0]
            sequence = lines[i + 1].strip()
            
            if (sequence[0:4] == 'CCTT' and sequence[18:22] == 'GAGC' and sequence[42:46] == 'ATAA'):
                capture_seq = 'CS2'
                barcode2 = sequence[22:42]
                read2_dict[identifier] = barcode2
                read2_cs_dict[identifier] = capture_seq

            elif (sequence[0:4] == 'TGCT' and sequence[17:21] == 'AAGC' and sequence[41:45] == 'ATAA'):
                capture_seq = 'CS1'
                barcode2 = sequence[21:41]
                read2_dict[identifier] = barcode2
                read2_cs_dict[identifier] = capture_seq
        
    # Write sequences to the output TSV file for matching identifiers
    with open(output_file, 'w') as file:
        file.write("umi\tbarcode1\tbarcode2\tread1_cs\tread2_cs\n")
        for identifier in set(read1_dict.keys()) & set(read2_dict.keys()) & set(read1_umi_dict.keys()) & set(read1_cs_dict.keys()) & set(read2_cs_dict.keys()):
            umi = read1_umi_dict[identifier]
            barcode1 = read1_dict[identifier]
            barcode2 = read2_dict[identifier]
            read1_cs = read1_cs_dict[identifier]
            read2_cs = read2_cs_dict[identifier]

            if (read1_cs == 'CS2' and read2_cs == 'CS1'):
                cs2_cs1_count+=1
            elif (read1_cs == 'CS1' and read2_cs == 'CS2'):
                cs1_cs2_count+=1
            elif (read1_cs == 'CS1' and read2_cs == 'CS1'):
                cs1_cs1_count+=1
            elif (read1_cs == 'CS2' and read2_cs == 'CS2'):
                cs2_cs2_count+=1

            file.write(f"{umi}\t{barcode1}\t{barcode2}\t{read1_cs}\t{read2_cs}\n")

    print(f"{name}\ncs2_cs1_count:{cs2_cs1_count}\ncs1_cs2_count:{cs1_cs2_count}\ncs1_cs1_count:{cs1_cs1_count}\ncs2_cs2_count:{cs2_cs2_count}\n\n\n",flush=True)

input_file ='samples.txt'
path = 'path/to/fastq'

with open(input_file,'r') as input_f:
	lines = input_f.readlines()

for i in lines:
		name = str(i.strip())
		read1_file = path+name+'_R1_001.fastq.gz'
		read2_file = path+name+'_R2_001.fastq.gz'																															
		output_file = name+'_extracted_bcs.tsv'

		extract_sequences(read1_file, read2_file, output_file)

