#$ -cwd
#$ -o run_cmp_cpu.output.txt
#$ -pe smp 1
#$ -q offline@cmpgrid-65
# NOTE: this script is executed by /bin/sh

date
. "/mnt/home.stud/lucivpav/.shrc"
cd "/mnt/datagrid/personal/lucivpav/InLocCIIRC_demo"
module load MATLAB/9.7
module load SuiteSparse/5.1.2-foss-2018b-METIS-5.1.0
module load LLVM/6.0.0-GCCcore-7.3.0 # avoid pyrender/OSMesa crash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/datagrid/personal/lucivpav/gflags-2.2.2/build/lib:/mnt/datagrid/personal/lucivpav/InLocCIIRC_demo/functions/vlfeat/toolbox/mex/mexa64
export PYTHONPATH=/mnt/home.stud/lucivpav/.local/lib/python3.5/site-packages
export INLOC_HW="CPU"
cat startup.m inloc_demo.m | matlab -nodesktop 2>&1 | tee InLocCIIRC_demo.log
date
