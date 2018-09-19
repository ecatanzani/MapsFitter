#!/bin/sh

run_exec () {
    TRY_IDX=$1
    BATCH_TRY=$2
    N_TRY=$3
    CMD="$EXEC ${TRY_IDX} ${BATCH_TRY} ${N_TRY}"; echo ${CMD}; ${CMD};
}

echo "`whoami` @ `hostname` : `pwd`"

CMD="source $SETVAR"; echo ${CMD}; ${CMD};

echo " ------------------------ "
echo " ---------- ENVS -------- "
echo " ------------------------ "
echo ""
env
echo ""
echo " ------------------------ "
echo " ------------------------ "
echo ""

HOME=`pwd`
UNIQUE="$$_$RANDOM"
SANDBOX=$HOME/../SandBox/$UNIQUE
while [[ -d $SANDBOX ]]; do
    UNIQUE="$$_$RANDOM"
    SANDBOX=$HOME/../SandBox/$UNIQUE
    echo "Searching for a free SandBox dir !"
done
mkdir -vp $SANDBOX
cd $SANDBOX
pwd

#Passing parameters to the function using a list file

#PARAMETERS="$PAR"
#cat <<< $PARAMETERS >> ./lista.dat
#run_exec lista.dat
#rm -fv lista.dat

#Passing each parameter separately to the function

run_exec "${PAR}"


ls -altrh
mv -v *.root *.txt $OUTDIR
rm -rfv "$SANDBOX"
