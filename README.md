# genome-shuffle-seq
Scripts, pipelines and sequences associated with genome shuffle seq described at: https://www.biorxiv.org/content/10.1101/2024.01.22.576756v1

This repository contains:
1. Pipeline and analysis scripts for analyzing amplicon-seq data
   - custom python scripts to extract barcode sequences from fastq files for libraries generated using either the 2 or 4 primer strategy as described in the methods
   - a sample input text file with the names of fastq files to be processed

2. Pipeline and analysis scripts for analyzing T7 IVT-seq data for mapping insertion sites in both mESCs and K562s
    - IVT_map_sort_pipeline.sh is used to run the pipeline for analyzing and definining position of integrants from IVT-seq data.
    - This pipeline runs custom python scripts that are included within the folder

3. Pipeline and analysis scripts for analyzing scRNAseq data after T7 IST, including de novo calling of clonotypes and assingment of cells to clonotypes
      
4. Sequences associated with novel plasmids described in the study
   
