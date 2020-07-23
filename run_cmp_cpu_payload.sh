#$ -cwd
#$ -pe smp 1
# NOTE: this script is executed by /bin/sh

date
. "/mnt/home.stud/lucivpav/.shrc"
module load MATLAB/9.7
module load CeresSolver/1.14.0-fosscuda-2018b
#module load SuiteSparse/5.1.2-foss-2018b-METIS-5.1.0 # TODO: verify GV part still works, then remove this line
module load LLVM/6.0.0-GCCcore-7.3.0 # avoid pyrender/OSMesa crash
# TODO: verify GV part still works, then remove the next line:
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/datagrid/personal/lucivpav/gflags-2.2.2/build/lib:"${PWD}"/functions/vlfeat/toolbox/mex/mexa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"${PWD}"/functions/vlfeat/toolbox/mex/mexa64
export PYTHONPATH=/mnt/home.stud/lucivpav/.local/lib/python3.5/site-packages
export INLOC_HW="CPU"
cat startup.m inloc_demo.m | matlab -nodesktop 2>&1
date
