#!/usr/bin/env python3
import sys
import re
import gzip
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import glob

#This script extracts, barcode1, barcode2, umi and strand from R1 of shuffle IVT sequencing data


def extract_sequences(sample, input_folder, output_folder):
    
    read1_file = glob.glob(input_folder + f"*{sample}*_R1_001.fastq.gz")[0]
    read2_file = glob.glob(input_folder + f"*{sample}*_R2_001.fastq.gz")[0]
    output_file = output_folder + "/" + sample + "_R1_bc.tsv"
    
    print(f'read1 file {read1_file}')
    print(f'read2 file {read2_file}')
    print(f'output file {output_file}')
    
    bc1_dict = {}
    read1_umi_dict = {}
    bc2_dict = {}
    bc1_dict_rc = {}
    bc2_dict_rc = {}
    strand_dict = {}
    recomb_dict = {}

    seq_dict = {}

    total_line = 0
    barcode_line = 0
                
    # Read and store sequences from Read 1 based on conditions
    with gzip.open(read1_file, 'rt') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            identifier = lines[i].split(' ')[0]
            sequence = lines[i + 1].strip()

            total_line += 1 

            if (sequence[28:32] == 'GAGC' and sequence[52:56] == 'ATAA' and sequence[10:14] == 'CCTT' and sequence[82:86] == 'TTAT'):
                strand = 'top-CS2'
                umi = sequence[:10]
                barcode1 = sequence[32:52]
                barcode1_revcomp = Seq(barcode1).reverse_complement()
                barcode2 = sequence[86:106]
                barcode2_revcomp = Seq(barcode2).reverse_complement()
                recombinase = 'loxP'

                strand_dict[identifier] = strand
                read1_umi_dict[identifier] = umi
                bc1_dict[identifier] = barcode1
                bc1_dict_rc[identifier] = barcode1_revcomp
                bc2_dict[identifier] = barcode2
                bc2_dict_rc[identifier] = barcode2_revcomp
                recomb_dict[identifier] = recombinase
            
            elif (sequence[28:32] == 'AAGC' and sequence[52:56] == 'ATAA' and sequence[10:14] == 'TTGC' and sequence[82:86] == 'TTAT'):
                strand = 'bottom-CS1'
                umi = sequence[:10]
                barcode1 = sequence[32:52]
                barcode1_revcomp = Seq(barcode1).reverse_complement()
                barcode2 = sequence[86:106]
                barcode2_revcomp = Seq(barcode2).reverse_complement()
                recombinase = 'loxP'

                strand_dict[identifier] = strand
                read1_umi_dict[identifier] = umi
                bc1_dict[identifier] = barcode1
                bc1_dict_rc[identifier] = barcode1_revcomp
                bc2_dict[identifier] = barcode2
                bc2_dict_rc[identifier] = barcode2_revcomp
                recomb_dict[identifier] = recombinase
            
            elif (sequence[28:32] == 'GAGC' and sequence[52:56] == 'GTGG' and sequence[10:14] == 'CCTT' and sequence[100:104] == 'CCCA'):
                strand = 'top-CS2'
                umi = sequence[:10]
                barcode1 = sequence[32:52]
                barcode1_revcomp = Seq(barcode1).reverse_complement()
                barcode2 = sequence[104:]
                barcode2_revcomp = Seq(barcode2).reverse_complement()
                recombinase = 'attP'

                strand_dict[identifier] = strand
                read1_umi_dict[identifier] = umi
                bc1_dict[identifier] = barcode1
                bc1_dict_rc[identifier] = barcode1_revcomp
                bc2_dict[identifier] = barcode2
                bc2_dict_rc[identifier] = barcode2_revcomp
                recomb_dict[identifier] = recombinase


            elif (sequence[28:32] == 'AAGC' and sequence[52:56] == 'TGGG' and sequence[10:14] == 'TTGC' and sequence[100:104] == 'CCAC'):
                strand = 'bottom-CS1'
                umi = sequence[:10]
                barcode1 = sequence[32:52]
                barcode1_revcomp = Seq(barcode1).reverse_complement()
                barcode2 = sequence[104:]
                barcode2_revcomp = Seq(barcode2).reverse_complement()
                recombinase = 'attP'

                strand_dict[identifier] = strand
                read1_umi_dict[identifier] = umi
                bc1_dict[identifier] = barcode1
                bc1_dict_rc[identifier] = barcode1_revcomp
                bc2_dict[identifier] = barcode2
                bc2_dict_rc[identifier] = barcode2_revcomp
                recomb_dict[identifier] = recombinase

            elif (sequence[28:32] == 'GAGC' and sequence[52:56] == 'GGCC' and sequence[10:14] == 'CCTT' and sequence[94:98] == 'CCGG'):
                strand = 'top-CS2'
                umi = sequence[:10]
                barcode1 = sequence[32:52]
                barcode1_revcomp = Seq(barcode1).reverse_complement()
                barcode2 = sequence[98:118]
                barcode2_revcomp = Seq(barcode2).reverse_complement()
                recombinase = 'attB'

                strand_dict[identifier] = strand
                read1_umi_dict[identifier] = umi
                bc1_dict[identifier] = barcode1
                bc1_dict_rc[identifier] = barcode1_revcomp
                bc2_dict[identifier] = barcode2
                bc2_dict_rc[identifier] = barcode2_revcomp
                recomb_dict[identifier] = recombinase


            elif (sequence[28:32] == 'AAGC' and sequence[52:56] == 'CCGG' and sequence[10:14] == 'TTGC' and sequence[94:98] == 'GGCC'):
                strand = 'bottom-CS1'
                umi = sequence[:10]
                barcode1 = sequence[32:52]
                barcode1_revcomp = Seq(barcode1).reverse_complement()
                barcode2 = sequence[98:118]
                barcode2_revcomp = Seq(barcode2).reverse_complement()
                recombinase = 'attB'

                strand_dict[identifier] = strand
                read1_umi_dict[identifier] = umi
                bc1_dict[identifier] = barcode1
                bc1_dict_rc[identifier] = barcode1_revcomp
                bc2_dict[identifier] = barcode2
                bc2_dict_rc[identifier] = barcode2_revcomp
                recomb_dict[identifier] = recombinase

        file.close()

    #open read2 fastq file and grab entire sequence to store with barcodes as separate column
    with gzip.open(read2_file, 'rt') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            identifier = lines[i].split(' ')[0]
            sequence = lines[i + 1].strip()
            seq_dict[identifier] = sequence

        file.close()

    # Write sequences to the output TSV file for matching identifiers
    with open(output_file, 'w') as file:
        file.write("identifier\tumi\tbarcode1\tbarcode1_rc\tbarcode2\tbarcode2_rc\tstrand\trecombinase_site\tseq\n")
        for identifier in set(bc1_dict.keys()) & set(bc2_dict.keys()) & set(read1_umi_dict.keys()) & set(strand_dict.keys()) & set(seq_dict.keys()):
            
            identifier=identifier
            umi = read1_umi_dict[identifier]
            barcode1 = bc1_dict[identifier]
            barcode1_revcomp = bc1_dict_rc[identifier]
            barcode2 = bc2_dict[identifier]
            barcode2_revcomp = bc2_dict_rc[identifier]
            strand = strand_dict[identifier]
            recombinase_site = recomb_dict[identifier]
            seq = seq_dict[identifier]

            file.write(f"{identifier}\t{umi}\t{barcode1}\t{barcode1_revcomp}\t{barcode2}\t{barcode2_revcomp}\t{strand}\t{recombinase_site}\t{seq}\n")
            barcode_line+= 1

        file.close()

    print("sample name: %s, total line: %f, lines with barcode: %f" %(sample, total_line, barcode_line))
    print(sample+" done")

if __name__ == "__main__":
    sample = sys.argv[1]
    input_folder = sys.argv[2]
    output_folder = sys.argv[3]
    
    extract_sequences(sample, input_folder, output_folder)    
        
