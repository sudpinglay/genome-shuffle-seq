#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -m ae
#$ -l mfree=100G,h_rt=6:00:00:00
#$ -pe serial 1
#$ -tc 2
#$ -t 1-2:1

# -M lalannej@uw.edu


module load modules modules-init modules-gs
module load python/3.12.1
module load pysam

#module load tbb/2020_U2 
module load bowtie2/2.5.3
module load samtools/1.14
module load seqtk/1.4

#module load gcc/8.1.0
module load R/4.3.2 


module load pear/0.9.11 


export PATH=/net/shendure/vol10/projects/Samuel/nobackup/10X/cellranger-6.0.1:$PATH

set -e

LOOKUP_FILE="/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/shuffleBC/lookup_shuffleBC_FB_20240905.txt"

if [[ ! -r "${LOOKUP_FILE}" ]]; then
    echo "Cannot find ${LOOKUP_FILE}" >&2
    exit 1
fi

LOOKUP="$(awk -v SGE_TASK_ID="${SGE_TASK_ID}" '$1 == SGE_TASK_ID {print $2}' < "${LOOKUP_FILE}")"
if [[ "x${LOOKUP}" = "x" ]]; then
    echo "Task ${SGE_TASK_ID} failed to lookup task ID" >&2
    exit 1
fi

date_str=$(date "+20%y%m%d")



dir_get_bc_out="/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/shuffleBC/get_BC_outs/"

# # # # # # # # # # # # # # # # # # # # 
# Calling cellRanger 
# # # # # # # # # # # # # # # # # # # # 

echo "running cellRanger"



path_fastq="/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/shuffleBC/"${LOOKUP}
date_str=$(date "+20%y%m%d")

cellranger count --id=${LOOKUP}"_CRv6_"${date_str} \
                 --transcriptome=/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/refdata-gex-GRCh38-2020-A \
                 --fastqs=${path_fastq} \
                 --localmem=64 \
                 --localcores=8



# # # # # # # # # # # # # # # # # # # # 
# get barcode (processing of cellRanger bam)
# # # # # # # # # # # # # # # # # # # # 

echo "get_barcode from cellRanger bam"

#source ~/virt_env/bin/activate
#module load modules modules-init modules-gs
#module load python/3.7.7
#module load pysam


input_bam_file=${LOOKUP}"_CRv6_"${date_str}"/outs/possorted_genome_bam.bam"
output_file1_name=${dir_get_bc_out}${LOOKUP}"_get_R2_2BC_"${date_str}".txt"

# parsing the cellRanger bam file
python /net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/shuffleBC/get_R2_2BC_T7shuffleAtt_20240905.py \
    --input_bam ${input_bam_file}\
    --output_file ${output_file1_name} \
    --chimeric_threshold 0.20 


# # # # # # # # # # # # # # # # # # # # 
# clean up UMI (Hamming distance correction)
# # # # # # # # # # # # # # # # # # # # 
echo "cleaning up UMIs"
# removing the G only mBCs (spurious and slow down UMI clean up)
#output_file1_name_no_G=${dir_get_bc_out}${LOOKUP}"_get_bc_v3_no_G_"${date_str}".txt"
#awk '$2!="GGGGGGGGGGGGGGGGGG" { print $0 }' ${output_file1_name} > ${output_file1_name_no_G}

# cleaning up the UMIs
#module load gcc/8.1.0
#module load R/3.5.1
output_file2_name=${dir_get_bc_out}${LOOKUP}"_get_R2_2BC_umi_cleaned"${date_str}".txt"
Rscript --vanilla /net/shendure/vol1/home/lalannej/vol10_projects_JB/seq015_mEB_crisprQTL_miniPilot_GFP_amplicons/nobackup/demux_fastq/clean_up_UMI_counts_v3_20220126.R \
	${output_file1_name} \
	${output_file2_name} \
	mBC



