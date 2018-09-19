#!/bin/sh

BATCH_TRY=1                     #Number of try for each job
N_TRY=2                         #Number of total try

STATUS=1

######################### Dependency paths ##########################

SEEDS_PATH="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/produceSeeds/seeds.txt"

SCALED_PATH="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/scalingSoftware/scaled_reference_Isotropic_histos.root"
DAMPE_ISO_MAP="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Salomon/results/FullHistos.root"


########### Final template paths

TEMPLATES_ALLSKY="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/computeTemplates/results/AllSkyTemplates.root"
TEMPLATES_DAMPE="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/computeTemplates/results/DAMPETemplates.root"

############################ Functions ##############################

function compute_free()
{
    FREE=$MAXPEND
    compute_jobs
    let FREE-=$JOBS
}

function scale_reference()
{
    REFERENCE_PATH="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Salomon/results/FullHistos.root"
    SCALE_EXE="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/scalingSoftware/SReference"

    CMD="$SCALE_EXE $REFERENCE_PATH $SCALED_PATH"
    echo ${CMD}
    ${CMD}
    STATUS=$?
}

function get_seeds()
{
    SEEDS_EXE="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/produceSeeds/PSeeds"

    CMD="$SEEDS_EXE $SEEDS_PATH $N_TRY $BATCH_TRY"
    echo ${CMD}
    ${CMD}
    STATUS=$?
}

function compute_templates()
{
    TEMPLATES_EXE="/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/Local/ExeSW/assets/computeTemplates/Release/JTemplates"

    CMD="$TEMPLATES_EXE $DAMPE_ISO_MAP $SCALED_PATH";
    echo ${CMD}
    ${CMD}
    STATUS=$?

}

function compute_jobs()
{
    JOBS=`qstat $QUEUE 2>&1 | egrep "WAIT|RUN" | wc -l`
}

function compute_maxpend()
{
    THRESHOLD_FILE=".THRESHOLD_PENDING"
    SHARPNESS=10
    echo "1000" > $THRESHOLD_FILE
    MAXPEND=1000
    MAXPENDL=$MAXPEND
    MAXPENDDEC=$MAXPEND
    let MAXPENDDEC/=$SHARPNESS
    let MAXPENDL-=MAXPENDDEC
}

####################################################################

compute_maxpend

SUBMITTED=0
compute_free
compute_jobs

########## SEED INDEX

TRY_IDX=0     #Starting line for the seeds.txt reading function

#####################


################# DEPENDENCY CHECK #################

if [[ ! -f ${SEEDS_PATH} ]]; then
    get_seeds
fi

if [[ ! -f ${SCALED_PATH} ]]; then
    scale_reference
fi

if [ ! -f ${TEMPLATES_ALLSKY} ] || [ ! -f ${TEMPLATES_DAMPE} ]; then
    compute_templates
fi

#####################################################

if [[ $STATUS -eq 1 ]]; then
    while [ $TRY_IDX -lt $N_TRY ]
    do
        if [[ $SUBMITTED -gt $FREE ]]; then
            while [ $JOBS -gt $MAXPENDL ]; do                                # MAXPENDL (MAXPEND-20%) to avoid to "release" the 'wait' just for one event
                echo "$JOBS job waiting..."
                sleep 10
                compute_jobs
                compute_maxpend
            done
            SUBMITTED=0
            compute_free
            echo "s=$SUBMITTED l=$FREE"
        fi
        let SUBMITTED+=1
        ./submit_one_job.sh $TRY_IDX $BATCH_TRY $N_TRY
        let TRY_IDX+=1
    done
fi

