#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -m ae
#$ -l mfree=100G,h_rt=6:00:00:00
#$ -pe serial 2
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



LOOKUP_FILE="/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/lookup_GEx_SP_shuffleBC_revision_20240903.txt"

if [[ ! -r "${LOOKUP_FILE}" ]]; then
    echo "Cannot find ${LOOKUP_FILE}" >&2
    exit 1
fi

LOOKUP="$(awk -v SGE_TASK_ID="${SGE_TASK_ID}" '$1 == SGE_TASK_ID {print $2}' < "${LOOKUP_FILE}")"

if [[ "x${LOOKUP}" = "x" ]]; then
    echo "Task ${SGE_TASK_ID} failed to lookup task ID" >&2
    exit 1
fi


#PARENT_DIR="/net/shendure/vol10/projects/JBL/seq055SP_T7_shuffle_10x/nobackup/GEx_fastqs/"

#DIR_LIST=("$PARENT_DIR"/*)
#FASTQ_DIR_TO_PROCESS=${DIR_LIST[$SGE_TASK_ID - 1]}

fastq_path_ori="/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/"
fastq_path=${fastq_path_ori}${LOOKUP}

date_str=$(date "+20%y%m%d")

cellranger count --id=${LOOKUP}"_CRv6_"${date_str} \
	--fastqs=$fastq_path_ori \
	--sample=${LOOKUP} \
	--transcriptome=/net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0 \
	--localmem=64 \
	--localcores=8

# mm10: 
# hg38: /net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/refdata-gex-GRCh38-2020-A \
	
