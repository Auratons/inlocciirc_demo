#!/bin/bash
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

rm -f run_cmp_cpu.output.txt
qsub -v INLOC_EXPERIMENT_NAME="${INLOC_EXPERIMENT_NAME}" -q offline@${HOSTNAME} -o /mnt/datagrid/personal/lucivpav/InLocCIIRC_dataset/logs/InLocCIIRC_demo_cpu-"${INLOC_EXPERIMENT_NAME}".log run_cmp_cpu_payload.sh