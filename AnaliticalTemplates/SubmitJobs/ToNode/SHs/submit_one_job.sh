#!/bin/sh

USER="$(whoami)"

#ENV="${HOME}/.bash_profile"
#if [[ ! -f $ENV ]]; then echo "$ENV does not exist" >&2; exit 1; fi

SUBSET_DIR="DAMPE_Sky"

############ SandBox for JobResults into local data dir ############

UNIQUE="$$_$RANDOM"
SANDBOX=/storage/gpfs_data/dampe/users/ecatanzani/Data/DAMPE/MapsFitter/CNAFJobs/tmpJobs/$UNIQUE
    while [[ -d $SANDBOX ]]; do
        UNIQUE="$$_$RANDOM"
        SANDBOX=/storage/gpfs_data/dampe/users/ecatanzani/Data/DAMPE/MapsFitter/CNAFJobs/tmpJobs/$UNIQUE
        echo "Searching for a free SandBox dir !"
    done
mkdir -vp $SANDBOX
cd $SANDBOX
pwd

######################

ENV_FILE=".bash_profile"
ENV_PATH="${HOME}/${ENV_FILE}"
if [[ ! -f ${ENV_PATH} ]]; then echo "$ENV_PATH does not exist" >&2; exit 1; fi

WORKDIR="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitter/AnaliticalTemplates/SubmitJobs/ToNode"
if [[ ! -d $WORKDIR ]]; then echo "$WORKDIR does not exist" >&2; exit 1; fi

EXEC="${WORKDIR}/ExeSW/Release/JMapsFit"
if [[ ! -f $EXEC ]]; then echo "$EXEC does not exist" >&2 ; exit 1; fi

SUBMITDIR="${WORKDIR}/SHs"
if [[ ! -d $SUBMITDIR ]]; then echo "$SUBMITDIR does not exist" >&2; exit 1; fi

TEMPLATE="${SUBMITDIR}/submit.job.template"
if [[ ! -f $TEMPLATE  ]]; then echo "$TEMPLATE does not exist" >&2; exit; fi

OUTDIR="${SANDBOX}/results/${SUBSET_DIR}"
if [[ ! -d $OUTDIR ]]; then mkdir -pv ${OUTDIR}
else echo "output dir will be $OUTDIR"; fi

QUEUE=dampe_usr

PAR=$@

NAME="MapsFit"

JOB_DIR="${SANDBOX}/out/jobs"
if [[ ! -d $JOB_DIR ]]; then mkdir -pv ${JOB_DIR}; fi

LOG_DIR="${SANDBOX}/out/log"
if [[ ! -d $LOG_DIR ]]; then mkdir -pv ${LOG_DIR}; fi

echo "Copying ENVIRONMENT variable to storage dierctory: "
SETVAR=${WORKDIR}/
cp -v ${ENV_PATH} ${SETVAR}
SETVAR+=$ENV_FILE

JOBNAME="job_${NAME}_${RANDOM}"
JOB=${JOB_DIR}/${JOBNAME}.job
ERRFILE=${LOG_DIR}/${JOBNAME}.err
LOGFILE=${LOG_DIR}/${JOBNAME}.log

rm -fv ${ERRFILE}
rm -fv ${LOGFILE}
cp -v ${TEMPLATE}                           ${JOB}
sed -i "s%_SETVAR_%${SETVAR}%g"             ${JOB}
sed -i "s%_EXEC_%${EXEC}%g"                 ${JOB}
sed -i "s%_PAR_%${PAR}%g"                   ${JOB}
sed -i "s%_OUTDIR_%${OUTDIR}%g"             ${JOB}
chmod 777 ${JOB}

echo " ----- $NAME"
CMD="bsub -n 1 -q ${QUEUE} -u ${USER} -J ${JOBNAME} -oo ${LOGFILE} -e ${ERRFILE} ${JOB}"
echo ${CMD}
${CMD}
