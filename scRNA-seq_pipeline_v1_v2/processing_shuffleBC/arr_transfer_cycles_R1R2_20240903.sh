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
module load seqtk/1.4


LOOKUP_FILE="/net/shendure/vol8/projects/JBL/seq075_SP_shuffleBC_revision_20240903/nobackup/shuffleBC/lookup_shuffleBC_libs_20240903.txt"

if [[ ! -r "${LOOKUP_FILE}" ]]; then
    echo "Cannot find ${LOOKUP_FILE}" >&2
    exit 1
fi

LOOKUP="$(awk -v SGE_TASK_ID="${SGE_TASK_ID}" '$1 == SGE_TASK_ID {print $2}' < "${LOOKUP_FILE}")"

if [[ "x${LOOKUP}" = "x" ]]; then
    echo "Task ${SGE_TASK_ID} failed to lookup task ID" >&2
    exit 1
fi


seqtk trimfq -e 47 ${LOOKUP}_R1_001.fastq.gz | \
	gzip > ${LOOKUP}_w_47bpR1R2transfer_R1_001.fastq.gz

seqtk trimfq -b 28 ${LOOKUP}_R1_001.fastq.gz | \
	join <(zcat ${LOOKUP}_R2_001.fastq.gz | nl) <(cat - | nl) | \
	awk -F ' ' '{if ($2 ~ /^@VH00979/ ) {print $4" "$5;} else {print $2$3;}}' | \
	gzip > ${LOOKUP}_w_47bpR1R2transfer_R2_001.fastq.gz