#!/usr/bin/env python3
import sys
import re
import pandas as pd
import operator
from collections import Counter
from cigar import Cigar

#this script parses alignments and assigns them to BL6 or CAST alleles

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def clean_up_variant_call(input):
    bedfile = open(input, "r")
    shortname = input.split(".")[0]
    temp = open(shortname + ".variants.txt", 'w')
    temp.write("chrom\tstart\tend\tstrand\tread\tpos\tcigar\tseq\tumi\tqscore\tsnp_chrom\tsnp_start\tsnp_end\tsnp\tsnp_type\tallele\n")
    
    total_count = 0
    no_snp_count = 0
    snv_count = 0
    indel_count = 0
    nei_count_snp = 0
    bl6_count_snp = 0 
    cast_count_snp = 0
    nei_count_indel = 0
    bl6_count_indel = 0 
    cast_count_indel = 0

    for line in bedfile:
        
        chrom,start,end,read = line.split("\t")[0:4]
        strand = line.split("\t")[5]
        cigar = line.split("\t")[7]
        seq,qscore,filename,snp_chrom,snp_start,snp_end,snp = line.split("\t")[11:18]
        umi = seq[0:8]

        total_count += 1


        if snp_start != '-1':
            snp_type = filename.split('.')[0].split('_')[-1]
            ref_seq = snp.split('/')[0].split('_')[1]
            snp_seq = snp.split('/')[1]
            pos = abs(int(start)-int(snp_start))
            
            if len(seq) > abs(pos):
                
                index = cigar.find('S')
                
                if index != -1 and index < len(cigar) - 1:
                    adj = int(cigar[:index])
                    pos = abs(pos)+adj

                x = Cigar(cigar)
                y = list(x.items())

                count = 0
                indel_positions = {}
                
                if 'SNV' in snp_type:
                    snv_count+=1

                    for i in y:
                        val = list(i)[0]
                        count = count + val
                        aln = list(i)[1]

                        if aln == 'M':
                            pass

                        if aln == 'I' and count < abs(pos):
                            pos=abs(pos)+val

                        if aln == 'D' and count < abs(pos):
                            pos=abs(pos)-val

                    if seq[abs(pos)] == snp_seq:
                        allele='cast'
                        cast_count_snp += 1

                    elif seq[abs(pos)] == ref_seq:
                        allele = 'bl6'
                        bl6_count_snp += 1

                    else:
                        allele = 'neither'
                        nei_count_snp = nei_count_snp + 1


                elif 'Indel' in snp_type:
                    indel_count+=1
                    foo = 0
                    for i in y:
                        val = list(i)[0]
                        count = count + val
                        aln = list(i)[1]
                        
                        if aln == 'M':
                            pass
                        
                        if aln == 'I':
                            if count < abs(pos):
                                pos=abs(pos)+val
                            
                            position = count-(val+1)
                            indel_positions[position] = val
                            
                        if aln == 'D':
                            if count < abs(pos):
                                pos=abs(pos)-val
                                foo+=1
                            
                            position = count-(val+1)
                            indel_positions[position] = val
                                                
                    indel_seq = seq[abs(pos):abs(pos)+len(snp_seq)]
                    indel_ref_seq = seq[abs(pos):abs(pos)+len(ref_seq)]
                    
                    if (pos+foo in indel_positions.keys()):
                        
                        if len(ref_seq) > len(snp_seq):
                            variant_type = 'deletion'
                            
                            if (indel_seq == snp_seq) and ('D' in cigar) and (len(ref_seq)-1 == indel_positions[pos+foo]):
                                allele='cast'
                                cast_count_indel+=1

                            elif indel_ref_seq == ref_seq:
                                allele = 'bl6'
                                bl6_count_indel+=1

                            else:
                                allele ='neither'
                                nei_count_indel+=1

                        elif len(ref_seq) < len(snp_seq):
                            variant_type = 'insertion'
                            
                            if (indel_seq == snp_seq) and ('I' in cigar) and (len(snp_seq)-1 == indel_positions[pos+foo]):
                                allele='cast' 
                                cast_count_indel+=1               

                            elif indel_ref_seq == ref_seq:
                                allele = 'bl6'
                                bl6_count_indel+=1

                            else:
                                allele ='neither'
                                nei_count_indel+=1
                                                
                    elif indel_ref_seq == ref_seq:
                        allele = 'bl6'
                        bl6_count_indel+=1
                    
                    else:
                        allele ='neither'
                        nei_count_indel+=1
                
                else:
                    print(f'failing 156 here\n{line}')
                    #allele ='neither'
                    #nei_count_indel+=1

            else:
                allele ='noVariant'
                snp_type = 'noSNP'
                no_snp_count+=1

        else:
            allele = 'noVariant'
            snp_type = 'noSNP'
            no_snp_count+=1

        if strand == "+":
            loc = end
        elif strand == "-":
            loc = start
            seq = reverse_complement(seq)

        temp.write(f"{chrom}\t{start}\t{end}\t{strand}\t{read}\t{loc}\t{cigar}\t{seq}\t{umi}\t{qscore}\t{snp_chrom}\t{snp_start}\t{snp_end}\t{snp}\t{snp_type}\t{allele}\n")

    print(f"{bedfile.name} variant calling done\n")
    print(f"total number of lines: {total_count}\n\nnumber of reads not overlapping variant: {no_snp_count}\number of lines in snvs: {snv_count}\nnumber of lines in indels: {indel_count}\n")
    print(f"number of snv lines assigned to bl6: {bl6_count_snp}\nnumber of snv lines assigned to cast: {cast_count_snp}\nnumber of snv lines assigned to neither: {nei_count_snp}")
    print(f"number of indel lines assigned to bl6: {bl6_count_indel}\nnumber of indel lines assigned to cast: {cast_count_indel}\nnumber of indel lines assigned to neither: {nei_count_indel}")
    
    missing = total_count - (no_snp_count + nei_count_indel + nei_count_snp + bl6_count_snp + bl6_count_indel + cast_count_snp + cast_count_indel) 
    print(f"\nmissing reads: {missing}")

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