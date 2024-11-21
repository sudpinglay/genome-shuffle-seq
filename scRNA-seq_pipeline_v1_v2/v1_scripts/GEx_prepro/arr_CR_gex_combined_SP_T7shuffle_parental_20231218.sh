#$ -q shendure-short.q
#$ -cwd
#$ -S /bin/bash
#$ -o /net/shendure/vol1/home/lalannej/sge_logs
#$ -e /net/shendure/vol1/home/lalannej/sge_logs
#$ -m ae
#$ -l mfree=100G,h_rt=6:00:00:00
#$ -pe serial 2
#$ -tc 1
#$ -t 1-1:1

# -M lalannej@uw.edu



module load modules modules-init modules-gs
module load pcre2/10.35
module load hdf5/1.10.1
module load R/4.0.0
module load gcc/8.2.0
export PATH="/net/shendure/vol1/home/sgregala/miniconda3/bin:$PATH"
module load python/3.7.7
export PATH=/net/shendure/vol10/projects/Samuel/nobackup/10X/cellranger-6.0.1:$PATH



#PARENT_DIR="/net/shendure/vol10/projects/JBL/seq055SP_T7_shuffle_10x/nobackup/GEx_fastqs/"

#DIR_LIST=("$PARENT_DIR"/*)
#FASTQ_DIR_TO_PROCESS=${DIR_LIST[$SGE_TASK_ID - 1]}
#fastq_path_ori="/net/shendure/vol10/projects/JBL/seq046_scQer_scv2_20230712/nobackup/GEx_demux/AACWKKMM5/outs/fastq_path/AACWKKMM5/"
#fastq_path=${fastq_path_ori}${LOOKUP}


date_str=$(date "+20%y%m%d")

cellranger count --id="sample_parental_combined_runs_GEx_CRv6_"${date_str} \
	--fastqs=/net/shendure/vol10/projects/JBL/seq055SP_T7_shuffle_10x/nobackup/GEx_fastqs/GEx_parental,/net/shendure/vol10/projects/JBL/seq055SP_T7_shuffle_10x/nobackup/GEx_fastqs/GEx_parental_reseq \
	--transcriptome=/net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0 \
	--localmem=64 \
	--localcores=8