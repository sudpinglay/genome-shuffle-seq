
import argparse
import pysam
import collections
import sys
import re
from Bio.Seq import Seq




if __name__ == '__main__':

    parser = argparse.ArgumentParser('Script to generate a file with counts of all mutation barcodes found within each cell in a CROP-seq or similar experiment.')
    parser.add_argument('--input_bam', '-i', help='Position sorted BAM (or list of bams) from 10X pipestance.')
    parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
    parser.add_argument('--whitelist', required=False, default=None, help='Optional mutation barcode whitelist.')
    parser.add_argument('--chimeric_threshold', type=float, default=0.2, help='Threshold for calling a UMI non-chimeric.')
    parser.add_argument('--force_correction', type=int, help='Force correction to a specified edit distance. Mismatches that can map to different sequences will be ignored and left uncorrected.')

    args = parser.parse_args()

    #search_seq=args.search_seq
    #barcode_start=args.seq_start
    #barcode_length=args.barcode_length
    input_bam=args.input_bam
    chimeric_threshold=args.chimeric_threshold
    output_file_name=args.output_file

    
    ############################################################
    # tallying up UMIs observed in cells for each mBC
    ############################################################
    
    mutation_barcodes = {}
    
    read_number=0
    read_mapped=0
    read_no_corrected_cBC_umi=0
    read_no_search_seq_match=0


    for read in pysam.Samfile(input_bam):
        read_number += 1

        if not read.is_unmapped:
            read_mapped += 1
            continue

        # get various part of the read in the bam file
        seq = read.seq.upper()
        tags = dict(read.tags)
        cell = tags.get('CB', None)
        umi = tags.get('UB', None)

        # skip read if no corrected UMI or cBC from cell ranger
        if not cell or not umi:
            read_no_corrected_cBC_umi += 1
            continue
        
        barcode = "None"
        
        ## search for the expe
        #match_CS1 = re.search(r'AAGC(.{20})ATAA', seq)
        #match_CS2 = re.search(r'GAGC(.{20})ATAA', seq)
        
        if re.findall(r'TGAGC(.{20})ATAAC',seq[15:49]) and (re.findall(r'GTTAT(.{20})',seq[69:])):
            capture_seq = 'CS2'
            barcode1 = re.findall(r'TGAGC(.{20})ATAAC',seq[15:49])[0] #this is directly realBC1
            barcode2 = re.findall(r'GTTAT(.{20})',seq[69:])[0] #barcode2 can be reverse complemented - realBC2
            barcode2_revcomp = Seq(barcode2).reverse_complement()
            
            barcode = barcode1+"_"+barcode2_revcomp+"_"+capture_seq

            #read2_barcode1_dict[identifier] = barcode1
            #read2_barcode2_dict[identifier] = barcode2_revcomp
            #read2_cs_dict[identifier] = capture_seq
                                
        elif (re.findall(r'AAAGC(.{20})ATAAC',seq[15:49])) and (re.findall(r'GTTAT(.{20})',seq[69:])):
            capture_seq = 'CS1'
            barcode2 = re.findall(r'AAAGC(.{20})ATAAC',seq[15:49])[0] #this is directly realBC2
            barcode1 = re.findall(r'GTTAT(.{20})',seq[69:])[0] #barcode1 can be reverse complememented -realBC1
            barcode1_revcomp = Seq(barcode1).reverse_complement()
            
            barcode = barcode1_revcomp+"_"+barcode2+"_"+capture_seq

        
        
        ##if match_CS1 != "None" & match_CS2=="None":
        #if match_CS1:
        #    barcode = match_CS1.group()
        #    position_BC = match_CS1.start()
        #    capture_seq = "CS1"
        ##elif match_CS2 != "None" & match_CS1=="None":
        #elif match_CS2:
        #    barcode = match_CS2.group()
        #    position_BC = match_CS2.start()
        #    capture_seq = "CS2"
        
       
        
        if barcode != "None":
            
            # not implementing error correction from on-list.
            corrected_barcode = barcode
            
            # list of all barcodes already listed w/ cBC=cell
            barcodes_in_cell = mutation_barcodes.get(cell, dict())

            # add the corresponding UMI to cBC
            if corrected_barcode not in barcodes_in_cell:
                barcodes_in_cell[corrected_barcode] = []
            
            barcodes_in_cell[corrected_barcode].append(umi)

            # update full dictionary with updated list of mBC/umis
            mutation_barcodes[cell] = barcodes_in_cell
            
        else:
            read_no_search_seq_match += 1
      
      
    print("read number:",read_number)
    print("mapped reads (not expected):",read_mapped)
    print("reads without corrected cBC+UMI:",read_no_corrected_cBC_umi)
    print("non mapped reads w/ corrected cBC+UMI:",read_number-read_mapped-read_no_corrected_cBC_umi)
    print("reads without exact match in search seq:", read_no_search_seq_match, read_no_search_seq_match/(read_number-read_mapped-read_no_corrected_cBC_umi))    
        
    ############################################################
    # generate output file (filter out chimeric UMIs)
    ############################################################
        
    output_file=open(output_file_name, 'w')
    original_stdout = sys.stdout
    output_file.write('\t'.join(['cBC','mBC','n_reads','n_UMI','n_reads_filtered','n_UMI_filtered','list_reads_per_UMIs','list_reads_per_UMIs_filtered'])+'\n')

    for cell in mutation_barcodes:

        # get the list of all UMI (each appearance corresponding to a read), irrespective of mBC, in a given cell (for a cBC)
        all_umis = []
        for mBC in mutation_barcodes[cell]:
            all_umis.extend(mutation_barcodes[cell][mBC])

        # tally up read counts for each UMI (again, irrespective of associated mBC)
        umi_counts_all_bc = collections.Counter(all_umis)    

        # loop through barcodes
        for mBC in mutation_barcodes[cell]:

            # counts of reads and UMIs for a given mBC, generate the list of read counts per UMI
            umi_counts_bc=collections.Counter(mutation_barcodes[cell][mBC])
            n_UMIs=len(umi_counts_bc.keys())
            n_reads=sum(umi_counts_bc.values())
            list_umi_counts=list(umi_counts_bc.items())
            str_list_umi_counts=[umi+"_"+str(count) for umi,count in list_umi_counts]

            # Get a list of UMIs that fall below the threshold for being likely non-chimeric.
            # the idea here: for a given cBC, each UMI should in principle be associated with a unique mBC. 
            # the proportion below is the fraction of reads from a UMI associated with the mBC of interest (over all other). If this fraction is too low
            chimeric_umis = set()
            for umi in umi_counts_bc:
                if float(umi_counts_bc[umi]) / umi_counts_all_bc[umi] < chimeric_threshold:
                    chimeric_umis.add(umi)

            # Filter any chimeric UMIs before outputing final stats
            filtered_umis = [umi for umi in mutation_barcodes[cell][mBC] if umi not in chimeric_umis]
            filtered_umi_counts_bc = collections.Counter(filtered_umis)
            n_UMIs_filtered=len(filtered_umi_counts_bc.keys())
            n_reads_filtered=sum(filtered_umi_counts_bc.values())
            list_umi_counts_filtered=list(filtered_umi_counts_bc.items())
            str_list_umi_counts_filtered=[umi+"_"+str(count) for umi,count in list_umi_counts_filtered]

            # output to file: only print if non 0 filtered umis>0. 
            if n_UMIs_filtered>0:
               output_file.write('\t'.join([cell,str(mBC),str(n_reads),str(n_UMIs),str(n_reads_filtered),str(n_UMIs_filtered)]))
               output_file.write('\t')
               sys.stdout = output_file
               print(*str_list_umi_counts,sep=",",end="\t")
               print(*str_list_umi_counts_filtered,sep=",")
               sys.stdout = original_stdout


        
        
        
        
