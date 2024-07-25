export SCRIPT_BASE=/camp/home/hungm/nemo-pipelines/
PIP_DIR=`pwd -L`
export PIP_DIR="${PIP_DIR##$SCRIPT_BASE}"
export now=$(date +'%Y%m%d')
export LOG_DIR=${PRJ}/logs/${PIP_DIR}/${SLURM_JOB_ID}

mkdir -p ${LOG_DIR}
rm ${SCRIPT_BASE}/${PIP_DIR}/slurm*
cp -r ${SCRIPT_BASE}/${PIP_DIR}/* ${LOG_DIR}/
exec > ${LOG_DIR}/${SLURM_JOB_ID}.log 2>&1

echo $(date)
