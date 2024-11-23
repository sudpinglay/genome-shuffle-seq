#!/bin/bash
#$ -N python_10x
#$ -l mfree=10G
#$ -cwd
#$ -pe serial 1-8
#$ -o python_stderr.txt
#$ -e python_stderr.txt
#$ -m abe
#$ -M pinglay@uw.edu

module load modules modules-init modules-gs
module load python/3.9.19_spec

rawdata="/net/shendure/vol10/projects/sud/nobackup/20240424_K562Shuffle/"
barcode_extract="amplicon_BCExtract_20240906.py"

#extract barcodes and UMIs
echo "Extracting barcodes..."
ls $rawdata*K562*R1_001.fastq.gz | cut -d "/" -f 9 | cut -d "R" -f 1 | parallel --gnu "python ${barcode_extract} {} $rawdata ."