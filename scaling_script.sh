#!/bin/bash

# NOTE:
# When running 'weak', make sure to specify config 
# related to specific system size.

# NOTE 2:
# For weak scaling, FIRST generate configs, then 
# run cfg_gen, io, flow.

FOLDER=$1
SYSTEM=$2
GAUGECONFIG=$3
DRYRUN=$4

set EMPTY=""
set GAUGECONFIGARGUMENT=""

if [[ "${FOLDER}" = "${EMPTY}" ]]
then
    echo "Missing folder of setup configurations."
    echo "Args: [folder] [system] [gaugeconfig]"
    exit 0
fi

if [[ "${SYSTEM}" = "${EMPTY}" ]]
then 
    echo "No system specified."
    exit 0
fi

if [[ "${DRYRUN}" = "dryrun" ]]
then
    DRYRUN="--dryrun"
else
    DRYRUN=""
fi

# For when we need to generate configurations 
# for each of the weak scaling tests.
if [[ "${GAUGECONFIG}" == "weak" ]]
then
    GAUGECONFIGARGUMENT=""
elif [[ "${GAUGECONFIG}" != "${EMPTY}" ]]
then 
    # echo "Missing configuration to pass to createJobs."
    GAUGECONFIGARGUMENT="-lcfg $GAUGECONFIG"
fi

ADDITIONAL_ARGS=""
if [[ $(basename $FOLDER) = "io" ]]
then
    ADDITIONAL_ARGS="-NCf 10"
    GAUGECONFIGARGUMENT="-lcfgr $GAUGECONFIG"
fi

# Loops over files in given folder
for filename in $FOLDER/*.py;
do
    echo ""
    echo "python createJobs.py $DRYRUN load $filename -s $SYSTEM $GAUGECONFIGARGUMENT --ignore_tasks_per_node $ADDITIONAL_ARGS"
    python2 createJobs.py $DRYRUN load $filename -s $SYSTEM $GAUGECONFIGARGUMENT --ignore_tasks_per_node $ADDITIONAL_ARGS
done
