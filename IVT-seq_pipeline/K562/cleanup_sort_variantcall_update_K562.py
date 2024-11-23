#!/usr/bin/env python3
import sys
import re
import pandas as pd
import operator
from collections import Counter
from cigar import Cigar

#this script parses alignments and assigns them to a position

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def clean_up_variant_call(input):
    bedfile = open(input, "r")
    shortname = input.split(".")[0]
    temp = open(shortname + ".variants.txt", 'w')
    temp.write("chrom\tstart\tend\tstrand\tread\tpos\tcigar\tseq\tumi\tqscore\n")
    
    for line in bedfile:
        
        chrom,start,end,read = line.split("\t")[0:4]
        strand = line.split("\t")[5]
        cigar = line.split("\t")[7]
        seq = line.split("\t")[11]
        qscore = line.split("\t")[12]

        if strand == "+":
            loc = end
        elif strand == "-":
            loc = start
            seq = reverse_complement(seq)

        umi = seq[0:8]

        temp.write(f"{chrom}\t{start}\t{end}\t{strand}\t{read}\t{loc}\t{cigar}\t{seq}\t{umi}\t{qscore}\n")

    print(f"{bedfile.name} bed file sorting done\n")
    
    temp.close()
    bedfile.close()

def sort_data(input):
    shortname = input.split(".")[0]
    filename = shortname + ".variants.txt"
    df = pd.read_table(filename)
    #df.columns = ["barcode", "chr", "pos", "start", "end", "strand", "cigar", "UMI", "seq"]
    df = df.sort_values(["chrom", "pos", "strand", "read", "umi"], ascending = (True, True, True, True, True))
    
    df.to_csv(shortname + ".variants.sorted.txt", sep = "\t", header = True, index = False)
    

if __name__ == "__main__":
    input = sys.argv[1]
   
    clean_up_variant_call(input)
    sort_data(input)