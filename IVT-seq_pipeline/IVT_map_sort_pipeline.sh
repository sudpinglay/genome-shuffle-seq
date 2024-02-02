#!/bin/bash
#$ -N T7_lib
#$ -l mfree=70G
#$ -cwd
#$ -pe serial 8
#$ -o out_2.txt
#$ -e out_2.txt
#$ -m abe
#$ -M pinglay@uw.edu


module load modules modules-init modules-gs
module load python/3.7.7
module load cutadapt/2.5
module load bwa/0.7.17
module load samtools/1.9
module load bedops/2.4.35
module load bedtools/2.29.2


rawdata="/net/shendure/vol10/projects/sud/20230511_NextSeq2000_ShuffleIVT/fastq/"
index="/net/shendure/vol10/projects/sud/index/bwa/mm10/mm10"
barcode_extract="IVT_extractBCs.py"
sam_filter="align_filter.py"
bed_filter="cleanup_sort_variantcall_update.py"
group_collapse="ivt_clustered_groupcollapse_iterable.py"

#extract barcodes and UMIs from read1
echo "Extracting barcodes ..."
ls $rawdata*IVT*R1_001.fastq.gz | cut -d"_" -f 7,8,9,10,11 | parallel --gnu "python ${barcode_extract} {} $rawdata ."

#Trim 3'ITR from reads and keep discard untrimmed reads (with a relatively high error rate, but )
echo "Adaptor trimming ..."
ls $rawdata*IVT*R2_001.fastq.gz | cut -d"_" -f 7,8,9,10,11 | parallel --gnu "cutadapt --cores=4 --discard-untrimmed -e 0.2 -m 10 -a CCCTAGAAAGATA -o {}_ITR.fastq.gz" $rawdata"985_Lipo_Unsorted_500_{}_R2_001.fastq.gz"

#Map trimmed reads to genome
echo "Alignment ..."
mkdir mapped/
ls *_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "bwa mem -Y $index {}.fastq.gz > mapped/{}.sam.temp"

#Sort sam files
echo "Sorting Sam files ..."
ls *_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "samtools view -b -F 4 mapped/{}.sam.temp | samtools sort -o mapped/{}.sort.temp.sam"

#Check end of reads. Must contain TTAA and the must align (>4M in CIGAR)
echo "Filtering Sam files ..."
ls *_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "python ${sam_filter} mapped/{}.sort.temp.sam"

echo "Running sam2bed ..."
ls *_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "sam2bed < mapped/{}.sam > mapped/{}.bed"

echo "Intersecting with known variants..."
ls mapped/*_ITR.bed | cut -d"." -f 1 | parallel --gnu "cut -f1-13 {}.bed | grep '^chr' | sort -k 1,1 -k2,2n > {}.noflag.bed"  #removing flag columns, removing all rows that don't start with "chr" and sorting like variant bed files
ls mapped/*_ITR.bed | cut -d"." -f 1 | parallel --gnu "bedtools intersect -a {}.noflag.bed -b /net/shendure/vol10/projects/sud/Cast_SNVs-Indels/Cast_SNVs.sorted.bed /net/shendure/vol10/projects/sud/Cast_SNVs-Indels/Cast_Indels.sorted.bed -loj -wa -wb -filenames -sorted > {}.overlappingVariants.txt"

echo "Assigning reads to variants..."
ls mapped/*overlappingVariants.txt | cut -d"." -f 1 | parallel --gnu "python ${bed_filter} {}.overlappingVariants.txt"

echo "Collapsing reads into groups..."
ls mapped/*.variants.sorted.txt | cut -d"." -f 1 | parallel --gnu "python ${group_collapse} {}.variants.sorted.txt"

rm mapped/*temp*

#original name: IVT_map_sort_updated_nonparental.sh

