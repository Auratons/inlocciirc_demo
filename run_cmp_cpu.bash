#!/bin/bash
echo -n "Are you on cmpgrid-65? [yes/no]: "
read PROCEED
if [[ "$PROCEED" != "yes" ]]; then
    exit 0
fi
echo "Detecting all jobs on current node:"
qstat -u "*"
echo -n "Is the queue empty (except for your own jobs)? [yes/no]: "
read PROCEED
if [[ "$PROCEED" != "yes" ]]; then
    exit 0
fi
rm -f run_cmp_cpu.output.txt
qsub run_cmp_cpu_payload.sh
