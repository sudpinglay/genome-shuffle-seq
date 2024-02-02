import gzip
import pandas as pd
import numpy as np

def extract_sequences(read1_file, read2_file, output_file):
    read1_dict = {}
    read1_umi_dict = {}
    read2_dict = {}

    # Read and store sequences from Read 1 based on conditions
    with gzip.open(read1_file, 'rt') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            identifier = lines[i].split(' ')[0]
            sequence = lines[i + 1].strip()
            if (sequence[28:32] == 'GAGC' and sequence[52:56] == 'ATAA'):
                barcode1 = sequence[32:52]
                read1_dict[identifier] = barcode1
            if (sequence[10:14] == 'CCTT'):
                umi = sequence[:10]
                read1_umi_dict[identifier] = umi
    
    # Read and store sequences from Read 2 based on conditions
    with gzip.open(read2_file, 'rt') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            identifier = lines[i].split(' ')[0]
            sequence = lines[i + 1].strip()
            if (sequence[17:21] == 'AAGC' and sequence[41:45] == 'ATAA'):
                barcode2 = sequence[21:41]
                read2_dict[identifier] = barcode2
        
    # Write sequences to the output TSV file for matching identifiers
    with open(output_file, 'w') as file:
        file.write("umi\tbarcode1\tbarcode2\n")
        for identifier in set(read1_dict.keys()) & set(read2_dict.keys()) & set(read1_umi_dict.keys()):
            umi = read1_umi_dict[identifier]
            barcode1 = read1_dict[identifier]
            barcode2 = read2_dict[identifier]
            file.write(f"{umi}\t{barcode1}\t{barcode2}\n")


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
		print(i+" done")


