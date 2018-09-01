#!/bin/bash

FOLDER=$1
SYSTEM=$2
GAUGECONFIG=$3

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

if [[ "${GAUGECONFIG}" != "${EMPTY}" ]]
then 
    # echo "Missing configuration to pass to createJobs."
    GAUGECONFIGARGUMENT="-lcfg $GAUGECONFIG"
fi

ADDITIONAL_ARGS=""
if [[ $(basename $FOLDER) = "io" ]]
then
    ADDITIONAL_ARGS="-NCf 10"
fi

# Loops over files in given folder
for filename in $FOLDER/*.py;
do
    echo ""
    echo "python createJobs.py --dryrun load $filename -s $SYSTEM $GAUGECONFIGARGUMENT --ignore_tasks_per_node $ADDITIONAL_ARGS"
    python2 createJobs.py --dryrun load $filename -s $SYSTEM $GAUGECONFIGARGUMENT --ignore_tasks_per_node $ADDITIONAL_ARGS
    echo "TEMP EXIT!";exit 0
done
