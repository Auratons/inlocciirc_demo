#!/bin/bash

if [ -z "$1" ]; then
    echo "Please provide an argument specifying what script to run. (without file extension)"
    echo "Suggestions: inloc_demo, PE_stability_check"
    exit 1
fi
export SCRIPT_NAME="$1"

if [ "${HOSTNAME}" != "cmpgrid-65" ] && [ "${HOSTNAME}" != "cmpgrid-69" ] && [ "${HOSTNAME}" != "cmpgrid-71" ]; then
    echo "You are not on the correct CPU node!"
    exit 0
fi

echo "Detecting all jobs on current node:"
qstat -u "*"
echo -n "Is the queue empty (except for your own jobs)? [yes/no]: "
read PROCEED
if [[ "$PROCEED" != "yes" ]]; then
    exit 0
fi

echo -n "Have you checked htop for running interactive jobs? [yes/no]: "
read PROCEED
if [[ "$PROCEED" != "yes" ]]; then
    exit 0
fi

if [ -z "${INLOC_EXPERIMENT_NAME}" ]; then
    export SUFFIX=""
else
    export SUFFIX="-${INLOC_EXPERIMENT_NAME}"
fi

rm -f run_cmp_cpu.output.txt
qsub -v INLOC_EXPERIMENT_NAME="${INLOC_EXPERIMENT_NAME}" -v SCRIPT_NAME="${SCRIPT_NAME}" -q offline@${HOSTNAME} -o /mnt/datagrid/personal/lucivpav/InLocCIIRC_dataset/logs/"${SCRIPT_NAME}"_cpu"${SUFFIX}".log run_cmp_cpu_payload.sh