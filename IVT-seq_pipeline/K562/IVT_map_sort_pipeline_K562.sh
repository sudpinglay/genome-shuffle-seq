#!/bin/bash
#$ -N T7_lib
#$ -l mfree=20G
#$ -cwd
#$ -pe serial 1-4
#$ -o out_K562.txt
#$ -e out_K562.txt


module load modules modules-init modules-gs
module load python/3.7.7
module load cutadapt/2.5
module load bwa/0.7.17
module load samtools/1.9
module load bedops/2.4.35
module load bedtools/2.29.2


rawdata="path to rawdata fastqs"
index="path to bwa hg38 index"
barcode_extract="IVT_extractBCs_updatedJul24.py"
sam_filter="align_filter.py"
bed_filter="cleanup_sort_variantcall_update_K562.py"
group_collapse="ivt_clustered_groupcollapse_iterable_K562.py"

#extract barcodes and UMIs from read1
echo "Extracting barcodes ..."
ls $rawdata*R1_001.fastq.gz | cut -d"/" -f 11 | cut -d"_" -f 1,2,3,4,5 | parallel --gnu --jobs 1 "python ${barcode_extract} {} $rawdata ."

#Trim 3'ITR from reads and keep discard untrimmed reads (with a relatively high error rate, but )
echo "Adaptor trimming ..."
ls $rawdata*R2_001.fastq.gz | cut -d"/" -f 11 | cut -d"_" -f 1,2,3,4,5 | parallel --gnu "cutadapt --cores=4 --discard-untrimmed -e 0.2 -m 10 -a CCCTAGAAAGATA -o {}_ITR.fastq.gz" $rawdata"{}_R2_001.fastq.gz"

#Map trimmed reads to genome
echo "Alignment ..."
mkdir mapped/
ls *K562*_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "bwa mem -Y $index {}.fastq.gz > mapped/{}.sam.temp"

#Sort sam files
echo "Sorting Sam files ..."
ls *K562*_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "samtools view -b -F 4 mapped/{}.sam.temp | samtools sort -o mapped/{}.sort.temp.sam"

#Check end of reads. Must contain TTAA and the must align (>4M in CIGAR)
echo "Filtering Sam files ..."
ls *K562*_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "python ${sam_filter} mapped/{}.sort.temp.sam"

echo "Running sam2bed ..."
ls *K562*_ITR.fastq.gz | cut -d"." -f 1 | parallel --gnu "sam2bed < mapped/{}.sam > mapped/{}.bed"

echo "Renaming chromosomes in bed file..."

substitutions="hg38_chrom_renaming.txt"

# Create a temporary sed script
sed_script=$(mktemp)

# Convert the substitutions file to a sed script
while read -r old_value new_value; do
    echo "s/$old_value/$new_value/g" >> "$sed_script"
done < "$substitutions"

# Process each file that ends with ITR.bed
for filename in mapped/*K562*ITR.bed; do
    # Apply the sed script to the file
    sed -i -f "$sed_script" "$filename"
done

# Remove the temporary sed script
rm "$sed_script"

echo "Filtering bed file..."
ls mapped/*K562*_ITR.bed | cut -d"." -f 1 | parallel --gnu "cut -f1-13 {}.bed | grep '^chr' | sort -k 1,1 -k2,2n > {}.noflag.bed"  #removing flag columns, removing all rows that don't start with "chr" and sorting like variant bed files
ls mapped/*K562*_ITR.bed | cut -d"." -f 1 | parallel --gnu "python ${bed_filter} {}.noflag.bed"

echo "Collapsing reads into groups..."
ls mapped/*K562*.variants.sorted.txt | cut -d"." -f 1 | parallel --gnu --halt soon,fail=1 --jobs 1 "python ${group_collapse} {}.variants.sorted.txt"

rm mapped/*temp*

#original name: IVT_map_sort_updated_nonparental.sh

