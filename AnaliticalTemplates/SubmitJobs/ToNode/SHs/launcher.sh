#!/bin/bash

function compute_free()
{
    FREE=$MAXPEND
    compute_jobs
    let FREE-=$JOBS
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

########## SEED INDEX

TRY_IDX=0     #Starting line for the seeds.txt reading function

#####################

compute_maxpend
compute_free
compute_jobs

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
