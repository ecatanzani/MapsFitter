#!/bin/sh

cat << "EOF"
_,.
,` -.)
'( _/'-\\-.
/,|`--._,-^|            ,
\_| |`-._/||          ,'|
|  `-, / |         /  /
|     || |        /  /
`r-._||/   __   /  /
__,-<_     )`-/  `./  /
'  \   `---'   \   /  /
|           |./  /
/           //  /
\_/' \         |/  /
|    |   _,^-'/  /
|    , ``  (\/  /_
\,.->._    \X-=/^
(  /   `-._//^`
`Y-.____(__}
|     {__)
()`

_____ ______  _______ ______      ______ ______         _______
(_____)  ___ \(_______)  ___ \    / _____)  ___ \   /\  (_______)
_  | |   | |_____  | |   | |  | /     | |   | | /  \  _____
| | | |   | |  ___) | |   | |  | |     | |   | |/ /\ \|  ___)
_| |_| |   | | |     | |   | |  | \_____| |   | | |__| | |
(_____)_|   |_|_|     |_|   |_|   \______)_|   |_|______|_|

_____     _                   _           _           _
(_____)   | |                 | |         (_)         (_)
_  ___ | | _      ___ _   _| | _  ____  _  ___  ___ _  ___  ____
| |/ _ \| || \    /___) | | | || \|    \| |/___)/___) |/ _ \|  _ \
___| | |_| | |_) )  |___ | |_| | |_) ) | | | |___ |___ | | |_| | | | |
(____/ \___/|____/   (___/ \____|____/|_|_|_|_(___/(___/|_|\___/|_| |_|

EOF



#####################################################################################
#####################################################################################

source config


function compute_free()
{
    FREE=$MAXPEND
    compute_jobs
    let FREE-=$JOBS
}

function scale_reference()
{

    O_STAT=$1
    MULT=$2

    COMPILE="make"
    CLEAN="make distclean"

    ${COMPILE}
    ${CLEAN}

    CMD="$SCALE_EXE $REFERENCE_PATH $SCALED_PATH $BINX_SCALING $BINY_SCALING $O_STAT $MULT"
    ${CMD}

    STATUS=$?
}

function get_seeds()
{
    COMPILE="make"
    CLEAN="make distclean"

    ${COMPILE}
    ${CLEAN}

    CMD="$SEEDS_EXE $SEEDS_PATH $N_TRY $BATCH_TRY"
    ${CMD}

    STATUS=$?
}

function compute_templates()
{
    COMPILE="make rebuild"
    CLEAN="make clean"

    ${COMPILE}
    ${CLEAN}

    CMD="$TEMPLATES_EXE $DAMPE_ISO_MAP $SCALED_PATH";
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
compute_free
compute_jobs

for stat in ${STATISTICS[*]}
do
    for idx in {1..2}
    do

        if [[ ! -f ${SEEDS_PATH} ]]; then
            get_seeds
        fi

        scale_reference ${array[$stat]} $idx
        compute_templates

        COMPILE="$"

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



    done
done

